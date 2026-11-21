from __future__ import annotations

import argparse
import io
import math
from pathlib import Path

import numpy as np
import pandas as pd


IMPACT_WEIGHTS = {
    "HIGH": 3.0,
    "MODERATE": 1.5,
    "LOW": 0.5,
    "MODIFIER": 0.2,
}

SAMPLE_VCF = """##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t10450\trs1\tA\tG\t99\tPASS\tGENE=TP53;IMPACT=HIGH
chr1\t10520\trs2\tC\tT\t85\tPASS\tGENE=TP53;IMPACT=MODERATE
chr2\t33100\trs3\tG\tA\t40\tLowQual\tGENE=BRCA1;IMPACT=HIGH
chr3\t55200\trs4\tT\tC\t95\tPASS\tGENE=GPX4;IMPACT=MODERATE
chr4\t88900\trs5\tA\tC\t90\tPASS\tGENE=ACSL4;IMPACT=LOW
"""

SAMPLE_DEG = """gene,log2FoldChange,padj
TP53,2.45,0.0012
GPX4,-1.85,0.0040
BRCA1,3.10,0.0005
ACSL4,0.40,0.4500
SLC7A11,1.90,0.0100
"""


def parse_info(info_str: str) -> dict[str, str]:
    data: dict[str, str] = {}
    for part in info_str.split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            data[k.strip()] = v.strip()
    return data


def load_vcf(vcf_text_or_path: str | None, require_pass: bool = True) -> pd.DataFrame:
    content = Path(vcf_text_or_path).read_text(encoding="utf-8") if vcf_text_or_path else SAMPLE_VCF
    rows: list[dict[str, object]] = []

    for line in content.splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) < 8:
            continue
        chrom, pos, var_id, ref, alt, qual, flt, info = parts[:8]
        if require_pass and flt not in ("PASS", "."):
            continue

        info_map = parse_info(info)
        gene = info_map.get("GENE")
        impact = info_map.get("IMPACT", "MODIFIER").upper()
        if not gene:
            continue

        rows.append(
            {
                "gene": gene.strip(),
                "impact": impact,
                "impact_weight": IMPACT_WEIGHTS.get(impact, 0.2),
            }
        )

    return pd.DataFrame(rows)


def load_deg(deg_path: str | None) -> pd.DataFrame:
    if deg_path:
        path = Path(deg_path)
        sep = "\t" if path.suffix.lower() in {".tsv", ".txt"} else ","
        df = pd.read_csv(path, sep=sep)
    else:
        df = pd.read_csv(io.StringIO(SAMPLE_DEG))

    df = df.copy()
    df["gene"] = df["gene"].astype(str).str.strip()
    df["log2FoldChange"] = pd.to_numeric(df["log2FoldChange"], errors="coerce")
    df["padj"] = pd.to_numeric(df["padj"], errors="coerce")
    df = df.dropna(subset=["gene", "log2FoldChange", "padj"])
    df["abs_lfc"] = df["log2FoldChange"].abs()
    df["neg_log10_padj"] = df["padj"].apply(lambda p: -math.log10(max(float(p), 1e-300)))
    return df


def prioritize(
    vcf_df: pd.DataFrame,
    deg_df: pd.DataFrame,
    padj_cutoff: float = 0.05,
    lfc_cutoff: float = 1.0,
) -> pd.DataFrame:
    if vcf_df.empty:
        return pd.DataFrame()

    var_summary = (
        vcf_df.groupby("gene", as_index=False)
        .agg(
            variant_count=("gene", "count"),
            max_impact_weight=("impact_weight", "max"),
            impacts=("impact", lambda s: ",".join(sorted(set(s)))),
        )
    )

    deg_df = deg_df.copy()
    deg_df["significant"] = (deg_df["padj"] <= padj_cutoff) & (deg_df["abs_lfc"] >= lfc_cutoff)

    merged = pd.merge(deg_df, var_summary, on="gene", how="inner")
    if merged.empty:
        return merged

    variant_score = merged["max_impact_weight"] * np.log1p(merged["variant_count"])
    expr_score = np.where(
        merged["significant"],
        merged["abs_lfc"] * np.minimum(merged["neg_log10_padj"], 10.0) / 5.0,
        0.2 * merged["abs_lfc"],
    )
    merged["priority_score"] = (variant_score + expr_score).round(4)

    cols = [
        "gene",
        "priority_score",
        "impacts",
        "variant_count",
        "log2FoldChange",
        "padj",
        "significant",
    ]
    return merged[cols].sort_values("priority_score", ascending=False).reset_index(drop=True)


def main() -> None:
    parser = argparse.ArgumentParser(description="Simple VCF + DEG candidate gene prioritizer.")
    parser.add_argument("--vcf", help="Optional VCF file path (uses built-in sample if omitted)")
    parser.add_argument("--deg", help="Optional DEG CSV/TSV file path (uses built-in sample if omitted)")
    parser.add_argument("--out", default="prioritized_genes.csv", help="Output CSV path")
    args = parser.parse_args()

    vcf_df = load_vcf(args.vcf)
    deg_df = load_deg(args.deg)
    result = prioritize(vcf_df, deg_df)

    print(result.to_string(index=False))
    result.to_csv(args.out, index=False)
    print(f"\nSaved to {args.out}")


if __name__ == "__main__":
    main()  
