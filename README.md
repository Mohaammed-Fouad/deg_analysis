# DEG Analysis

A command-line pipeline that identifies **differentially expressed genes (DEGs)** between a diseased group and a healthy group from a gene-expression matrix, originally built around a **lung cancer** dataset, but not limited to it: any two-group expression matrix with matching sample metadata works.

The pipeline computes fold change and Welch's t-test statistics per gene (no DEG-calling package used, only NumPy/SciPy), applies Benjamini-Hochberg multiple-testing correction implemented from scratch, classifies every gene, draws a volcano plot, runs a co-expression analysis between up- and down-regulated genes, and runs GO Biological Process enrichment via Enrichr, all from a single non-interactive command.

> **Note:** GO enrichment (Enrichr) matches genes by *gene symbol*, not by *probe/feature ID*. If your expression matrix is indexed by probe IDs (e.g. Affymetrix `1007_s_at`), use `--id-map` to translate to gene symbols first, see [A note on gene identifiers and GO enrichment](#a-note-on-gene-identifiers-and-go-enrichment) below.

---

## Table of Contents

- [Overview](#overview)
- [How It Works](#how-it-works)
- [Requirements](#requirements)
- [Installation](#installation)
- [Usage](#usage)
  - [A note on gene identifiers and GO enrichment](#a-note-on-gene-identifiers-and-go-enrichment)
  - [All available options](#all-available-options)
- [Project Structure](#project-structure)
- [Output Files](#output-files)
- [Example Run](#example-run)
- [Notes and Limitations](#notes-and-limitations)
- [Acknowledgements](#acknowledgements)
- [License](#license)

---

## Overview

Given a gene-expression matrix and a sample metadata file that labels each sample as diseased or healthy, the pipeline:

1. **Loads and splits** the expression matrix into disease / healthy groups. Optionally maps probe/feature IDs to gene symbols via `--id-map`, since GO enrichment only matches on gene symbols.
2. **Computes** log2 fold change, Welch's t-test statistic, and p-value for every gene, fully vectorized across all genes at once.
3. **Adjusts p-values** for multiple testing with the Benjamini-Hochberg procedure, implemented from scratch with NumPy.
4. **Classifies** every gene as over-expressed, down-expressed, or equally-expressed.
5. **Draws** a volcano plot (fold change vs. -log10 adjusted p-value).
6. **Saves** a tab-separated DEG report (gene, fold_change, p_value, adj_p_value, t_stat, change_type).
7. **Runs a co-expression analysis** between every up-regulated gene and every down-regulated gene (Pearson correlation on the diseased samples), reporting statistically significant similar (positive) and different (negative) pairs, ranked by correlation.
8. **Runs GO Biological Process enrichment** (via `gseapy`/Enrichr) separately for the up-regulated genes, the down-regulated genes, and all DEGs together.

This is a straight command-line port of the original Jupyter-notebook prototype: the statistics, the BH correction, the classification rules, the volcano-plot mechanics, and the co-expression logic are all unchanged. What changed is packaging: everything now lives in importable functions driven by a proper CLI (`argparse`) with a `main()` entry point, so the whole pipeline can be run non-interactively, with logging and configurable parameters, instead of by hand cell-by-cell.

---

## How It Works

| Stage | Description | Key Library |
|---|---|---|
| Data loading | Load metadata + expression matrix, split into disease/healthy groups, optionally map probe IDs to gene symbols | `pandas` |
| Fold change + t-test | Compute log2 fold change and Welch's t-test statistic/p-value for every gene at once (vectorized) | `numpy`, `scipy` |
| Multiple-testing correction | Benjamini-Hochberg FDR adjustment, implemented from scratch | `numpy` |
| Classification | Label each gene over-/down-/equally-expressed by fold-change and adjusted-p cutoffs | `numpy` |
| Volcano plot | Fold change vs. -log10 adjusted p-value, colored by classification | `matplotlib` |
| DEG report | Tab-separated report of every gene's stats and classification | `pandas` |
| Co-expression analysis | Pearson correlation between every up-gene/down-gene pair, computed as one correlation matrix | `numpy`, `scipy` |
| GO enrichment | Enrichr GO Biological Process enrichment for up, down, and all DEGs | `gseapy` |

---

## Requirements

- **Python 3.9+**
- An internet connection (for the GO enrichment step, unless `--skip-go` is used).

Install everything with a single command:

```bash
pip install -r requirements.txt
```

`requirements.txt`:

```
numpy>=1.23
pandas>=1.5
scipy>=1.9
matplotlib>=3.6
gseapy>=1.0
```

---

## Installation

1. **Clone the repository:**

   ```bash
   git clone https://github.com/Mohaammed-Fouad/deg-analysis.git
   cd deg-analysis
   ```

2. **(Recommended) Create a virtual environment:**

   ```bash
   python3 -m venv venv
   source venv/bin/activate   # Linux/Mac
   venv\Scripts\activate      # Windows PowerShell
   ```

3. **Install all dependencies:**

   ```bash
   pip install -r requirements.txt
   ```

---

## Usage

### 1. Run the pipeline with default settings

```bash
python deg_analysis.py --metadata meta.csv --expression expr.csv --outdir results
```

`--metadata` and `--expression` are always required. By default, samples are split using `group` column values `tumor` (disease) and `normal` (healthy).

### 2. Customize group labels and DEG thresholds

```bash
python deg_analysis.py --metadata meta.csv --expression expr.csv --disease-label tumor --healthy-label normal --fc-cutoff 1.0 --pval-cutoff 0.05 --skip-go
```

`--skip-go` is useful when running offline or when the GO enrichment step isn't needed.

### 3. Map probe IDs to gene symbols before GO enrichment

```bash
python deg_analysis.py --metadata meta.csv --expression expr.csv --id-map GPL570-55999.txt --id-map-id-col ID --id-map-symbol-col "Gene Symbol"
```

### A note on gene identifiers and GO enrichment

Enrichr (used for the GO enrichment step) matches submitted genes by *gene symbol* (e.g. `TP53`, `EGFR`). Many microarray expression matrices are indexed by *probe/feature ID* instead (e.g. Affymetrix IDs like `1007_s_at`), which will never match anything in Enrichr's gene sets. **This fails silently:** the API call still "succeeds," it just returns an empty results table (0 annotations) for every gene list.

If your gene index looks like probe IDs, use `--id-map` to supply a CSV that maps each probe/feature ID to its gene symbol, e.g. a GEO platform annotation file, which for many Affymetrix platforms already uses `ID` and `Gene Symbol` as column names, the defaults for `--id-map-id-col` / `--id-map-symbol-col`. The script also logs a warning if it detects Affymetrix-style probe IDs and no `--id-map` was given.

### All available options

```bash
python deg_analysis.py --help
```

**Input files**

| Flag | Description | Default |
|---|---|---|
| `--metadata` | CSV file with sample metadata (sample IDs + group labels) | **required** |
| `--expression` | CSV file with the gene expression matrix (genes as rows, samples as columns; first column = gene name) | **required** |
| `--sample-col` | Column in the metadata file holding sample IDs | `sample_accession` |
| `--group-col` | Column in the metadata file holding the condition label | `group` |
| `--disease-label` | Value in `--group-col` that marks a diseased sample | `tumor` |
| `--healthy-label` | Value in `--group-col` that marks a healthy/control sample | `normal` |
| `--id-map` | Optional CSV mapping probe/feature IDs to gene symbols (e.g. a GEO platform annotation file) | *(none)* |
| `--id-map-id-col` | Column in `--id-map` holding the original probe/feature ID | `ID` |
| `--id-map-symbol-col` | Column in `--id-map` holding the gene symbol to map to | `Gene Symbol` |
| `--id-map-agg` | How to aggregate multiple probes mapping to the same gene symbol: `mean`, `median`, or `max` | `mean` |

**DEG calling thresholds**

| Flag | Description | Default |
|---|---|---|
| `--fc-cutoff` | Absolute log2 fold-change cutoff for calling a gene a DEG | `1.0` |
| `--pval-cutoff` | Adjusted p-value cutoff for calling a gene a DEG | `0.05` |
| `--corr-pval-cutoff` | P-value cutoff for calling a co-expression pair significant | `0.05` |

**GO enrichment**

| Flag | Description | Default |
|---|---|---|
| `--skip-go` | Skip the GO enrichment step (e.g. no internet access) | off |
| `--gene-sets` | `gseapy`/Enrichr gene-set library to use | `GO_Biological_Process_2021` |
| `--organism` | Organism passed to `gseapy.enrichr` | `human` |

**Output**

| Flag | Description | Default |
|---|---|---|
| `--outdir` | Directory where all output files are written | `deg_output` |
| `--dump-all-frames` | Also dump every intermediate DataFrame the pipeline computes (not just the curated reports) as CSVs under `<outdir>/all_dataframes/` | off |
| `-v`, `--verbose` | Enable verbose (DEBUG) logging | off |

---

## Project Structure

```
.
├── deg_analysis.py     # Main pipeline script (single entry point)
├── requirements.txt    # One-command dependency installation
├── README.md           # This file
└── <outdir>/            # Generated: volcano plot, reports, GO enrichment results
```

---

## Output Files

Running the pipeline generates, under `--outdir`:

- `volcano_plot.png`: fold change vs. -log10 adjusted p-value, colored by classification.
- `deg_report.txt`: tab-separated report with one row per gene (`gene`, `fold_change`, `p_value`, `adj_p_value`, `t_stat`, `change_type`).
- `coexpression_similar.tsv` / `coexpression_different.tsv`: significant positively-/negatively-correlated up-gene/down-gene pairs, ranked by correlation strength.
- `go_enrichment_up/`, `go_enrichment_down/`, `go_enrichment_all/`: `gseapy`'s per-list Enrichr output folders (skipped entirely if `--skip-go` is set).

If `--dump-all-frames` is set, `<outdir>/all_dataframes/` additionally contains: the full post-`--id-map` expression matrix, the disease and healthy sub-matrices, the full-precision (unrounded) DEG stats table, the full-precision co-expression pair tables, and each Enrichr results table (if GO enrichment ran).

---

## Example Run

```
(venv) PS D:\deg_analysis> python deg_analysis.py --metadata meta.csv --expression expr.csv --id-map GPL570-55999.txt --id-map-id-col ID --id-map-symbol-col "Gene Symbol" --dump-all-frames --outdir results
01:55:35 [INFO] Loading metadata from meta.csv
01:55:35 [INFO] Found 60 disease samples and 60 healthy samples.
01:55:35 [INFO] Loading expression matrix from expr.csv
01:55:35 [INFO] Loading gene ID mapping from GPL570-55999.txt
01:55:36 [INFO] Loaded 45782 probe/feature ID -> gene symbol mappings.
01:55:36 [WARNING] Dropped 8893/54675 rows with no gene symbol in --id-map.
01:55:36 [INFO] Collapsed 45782 probes mapping to duplicate gene symbols into 23520 unique genes using 'mean' aggregation.
01:55:36 [INFO] Diseased matrix shape: (23520, 60) | Healthy matrix shape: (23520, 60)
01:55:39 [INFO] Detected that the expression matrix is ALREADY log-transformed.
01:55:40 [INFO] Gene regulation breakdown:
change_type
equally-expressed    22654
down-expressed         577
over-expressed         289
01:55:40 [INFO] Volcano plot saved to results\volcano_plot.png
01:55:40 [INFO] DEG report saved to results\deg_report.txt
01:55:41 [INFO] Co-expression reports saved to results\coexpression_similar.tsv and results\coexpression_different.tsv
01:55:41 [INFO] Found 5106 significant similar (positive) pairs and 54331 significant different (negative) pairs.
01:55:41 [INFO] Sending profiles to Enrichr API: 289 up, 577 down genes...
01:55:43 [INFO] Success: up pathway processing saved. Found 1607 annotations.
01:55:48 [INFO] Success: down pathway processing saved. Found 2893 annotations.
01:55:54 [INFO] Success: all pathway processing saved. Found 3438 annotations.
01:55:57 [INFO] Done. All outputs written to D:\deg_analysis\results
01:55:57 [INFO] All intermediate dataframes dumped to D:\deg_analysis\results\all_dataframes
```

---

## Notes and Limitations

- No dedicated DEG-calling package is used: fold change, Welch's t-test, and the Benjamini-Hochberg correction are all implemented from scratch with NumPy/SciPy.
- The pipeline auto-detects whether the expression matrix is already log-transformed by checking the maximum value (below 25 is treated as already logged); double-check the logged `[INFO]` line matches your data if this matters for your analysis.
- Co-expression p-values are derived analytically from the Pearson correlation's t-distribution, not via permutation, so results assume the standard parametric assumptions hold.
- GO enrichment requires an internet connection and depends on the Enrichr API's availability; use `--skip-go` to run fully offline.
- If your gene index is probe/feature IDs rather than gene symbols, GO enrichment will silently return 0 annotations unless you supply `--id-map`, see [A note on gene identifiers and GO enrichment](#a-note-on-gene-identifiers-and-go-enrichment).
- This tool is intended for research and educational purposes and does **not** constitute medical or diagnostic advice.

---

## Acknowledgements
This pipeline was originally developed as a term project for the **CIT672: Programming for Bioinformatics** course at **Nile University**.
### Group (A) Members
1. [Akram Nader Salah](https://github.com/akramnader289)
2. [Mariam Mohamed Mahdy](https://github.com/MariamMahdy)
3. [Ala'a Ahmed Fathallah](https://github.com/AlaaAzzam2311)
4. [Mohammed Fouad Mohammed](https://github.com/Mohaammed-Fouad)

**Project Topic:** *Differentially Expressed Genes (DEG) Analysis for Lung Cancer and Gene Ontology (GO) Enrichment*

> Thank you to everyone in the group for their collaboration, dedication, and shared efforts throughout this project.
> A special thanks to our course instructor, **Dr. Ibrahim Mohamed Youssef**, for his guidance, effort, and support during the course.

---

## License

No license has been specified for this project yet. If you plan to share or reuse this code, add a `LICENSE` file (e.g. [MIT](https://opensource.org/licenses/MIT)) to make the terms explicit.
