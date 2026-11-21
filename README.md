# var-deg-prioritizer
Multi-omics Python CLI tool that merges filtered VCF variants and RNA-seq differential expression data to rank and visualize candidate genes.
# var-deg-prioritizer



`var-deg-prioritizer` combines variant-level evidence from VCF files (`FILTER`, variant impact annotations such as `HIGH` or `MODERATE`, and per-gene variant counts) with gene-level differential expression statistics (`log2FoldChange` and `padj`) to compute a transparent composite `priority_score`.

## Project Structure

```text
var-deg-prioritizer/
├── README.md
├── requirements.txt
├── data/
│   └── sample/
│       ├── variants.vcf
│       └── deg_table.csv
├── src/
│   └── var_deg_prioritizer/
│       ├── __init__.py
│       ├── io_utils.py
│       ├── scorer.py
│       ├── plotter.py
│       └── cli.py
└── tests/
    └── test_scorer.py
