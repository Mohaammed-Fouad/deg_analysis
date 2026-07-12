# DEG Analysis CLI

A command-line pipeline that identifies **Differentially Expressed Genes (DEGs)** between a diseased sample group and a healthy/control sample group from a gene-expression matrix, then characterizes those DEGs through co-expression analysis and Gene Ontology (GO) enrichment.

The pipeline takes a gene-expression matrix and a sample metadata file as input (e.g. a tumor-vs-normal microarray or RNA-seq dataset such as GSE19804), computes fold change and statistical significance for every gene from scratch with NumPy/SciPy (no package infers the DEGs for you), corrects for multiple testing, draws a volcano plot, and reports up- and down-regulated genes together with their co-expression partners and enriched biological functions.

> **Note:** there is no default disease/control dataset bundled with this tool. You must always point it at your own `--metadata` and `--expression` files.

---

## Table of Contents

- [Overview](#overview)
- [How It Works](#how-it-works)
- [Requirements](#requirements)
- [Installation](#installation)
- [Input File Formats](#input-file-formats)
- [Usage](#usage)
- [Statistical Methodology](#statistical-methodology)
- [Project Structure](#project-structure)
- [Output Files](#output-files)
- [Troubleshooting](#troubleshooting)
- [Notes and Limitations](#notes-and-limitations)
- [License](#license)

---

## Overview

Given a gene-expression matrix and a metadata file labeling each sample as diseased or healthy, the pipeline:

1. **Loads** the expression matrix and metadata, and splits the matrix into a disease-group sub-matrix and a healthy-group sub-matrix based on the labels you specify.
2. **Computes per-gene statistics** - log2 fold change and a Welch's t-test (statistic + p-value) for every gene, fully vectorized across the whole matrix at once rather than looping gene by gene.
3. **Adjusts p-values** for multiple testing using a from-scratch, vectorized implementation of the Benjamini-Hochberg (FDR) procedure - no external multiple-testing package is used.
4. **Classifies** every gene as `over-expressed`, `down-expressed`, or `equally-expressed` against configurable fold-change and adjusted-p-value cutoffs.
5. **Draws a volcano plot** (fold change vs. -log10 adjusted p-value) with matplotlib.
6. **Saves a DEG report** - one line per gene: name, fold change, p-value, adjusted p-value, t-statistic, and change type.
7. **Runs a co-expression analysis** between every up-regulated and down-regulated gene (Pearson correlation on the diseased samples), reporting statistically significant positive ("similar") and negative ("different") pairs, ranked by correlation strength.
8. **Runs GO Biological Process enrichment** (via `gseapy`/Enrichr) separately for the up-regulated genes, the down-regulated genes, and all DEGs together, so you can compare which biological functions each subset is enriched for.

---

## How It Works

| Stage | Description | Key Library |
|---|---|---|
| Data loading | Read metadata + expression matrix, split into disease/healthy sub-matrices by sample group | `pandas` |
| Per-gene statistics | Vectorized log2 fold change and Welch's t-test across every gene at once | `numpy`, `scipy.stats` |
| Multiple-testing correction | Vectorized Benjamini-Hochberg FDR adjustment (`argsort` / `minimum.accumulate` / `clip`) | `numpy` |
| Classification | Label genes over-/down-/equally-expressed via cutoffs (`np.select`) | `numpy` |
| Volcano plot | Fold change vs. -log10 adjusted p-value, color-coded by change type | `matplotlib` |
| DEG report | Tab-separated report of every gene's stats and call | `pandas` |
| Co-expression analysis | Up-vs-down gene correlation matrix (`np.corrcoef`) with analytic p-values from the t-distribution | `numpy`, `scipy.stats` |
| GO enrichment | Enrichr GO Biological Process enrichment for up-only, down-only, and all DEGs | `gseapy` |

---

## Requirements

- **Python 3.9+**
- An internet connection for the GO enrichment step (it calls the public Enrichr API) - pass `--skip-go` to run fully offline.

Everything else - data loading, statistics, plotting, co-expression - runs locally with no external services.

---

## Installation

1. **Clone the repository** (or just download `deg_analysis.py`):

   ```bash
   git clone https://github.com/Mohaammed-Fouad/deg_analysis.git
   cd deg_analysis
   ```

2. **(Recommended) Create a virtual environment:**

   ```bash
   python3 -m venv venv
   source venv/bin/activate   # Linux/Mac
   venv\Scripts\activate      # Windows PowerShell
   ```

3. **Install all dependencies with a single command:**

   ```bash
   pip install -r requirements.txt
   ```

---

## Input File Formats

**Expression matrix** (`--expression`): CSV with genes as rows and samples as columns. The first column must be the gene name/ID (used as the row index).

| gene | T1 | T2 | N1 | N2 |
|--------|------|------|------|------|
| GENE1 | 8.21 | 7.98 | 5.10 | 5.32 |
| GENE2 | 4.55 | 4.61 | 4.50 | 4.48 |

**Sample metadata** (`--metadata`): CSV mapping each sample ID to a condition/group label.

| sample_accession | group |
|-------------------|--------|
| T1 | tumor |
| T2 | tumor |
| N1 | normal |
| N2 | normal |

Column names and label values don't have to match this example - point the script at whatever columns/labels your files actually use with `--sample-col`, `--group-col`, `--disease-label`, and `--healthy-label`.

### A note on gene identifiers and GO enrichment

GO enrichment (`gseapy`/Enrichr) matches submitted genes by **gene symbol** (e.g. `TP53`, `EGFR`). Many microarray platforms - including the Affymetrix arrays commonly used for datasets like GSE19804 - index their expression matrix by **probe/feature ID** instead (e.g. `1007_s_at`), not gene symbol. Submitting probe IDs to Enrichr doesn't raise an error; it just silently returns **0 annotations** for every gene list, because none of the IDs match anything in Enrichr's gene sets.

The script detects this automatically: if you don't pass `--id-map` and the gene index looks like Affymetrix-style probe IDs, it logs a warning up front. If a GO enrichment call still comes back with 0 annotations, it's also flagged explicitly in the log rather than reported as a plain "Success."

**To fix it, map probes to gene symbols before running the pipeline** with `--id-map`:

```bash
python deg_analysis.py --metadata meta.csv --expression expr.csv --id-map GPL570-55999.txt --id-map-id-col ID --id-map-symbol-col "Gene Symbol"
```

`--id-map` accepts any CSV/TSV with an ID column and a gene-symbol column - a GEO platform annotation file (e.g. `GPLxxxx-yyyyy.txt`, downloadable from the platform's page on [GEO](https://www.ncbi.nlm.nih.gov/geo/)) works directly, since many Affymetrix platform files already use `ID` and `Gene Symbol` as column headers (the tool's defaults). Probes with no symbol are dropped, and multiple probes mapping to the same gene symbol are collapsed into one row via `--id-map-agg` (`mean` by default; `median`/`max` also available).

### Getting and cleaning a GPL annotation file (common pitfalls)

Datasets like GSE19804 don't ship a ready-made ID-to-symbol file - you have to pull one from GEO yourself, and the file GEO gives you usually isn't clean enough to load right away. Two problems come up in practice:

**1. You don't already have an `--id-map` file - it isn't part of the expression matrix.**
It's tempting to point `--id-map` at your `*_expression_matrix.csv`, since it also has an "ID" column (`ID_REF`) - but that file only has probe IDs and per-sample expression values, never gene symbols. `--id-map` must be a *separate* file: the annotation table for the microarray **platform** the samples were run on, not the samples themselves.

To find it:
- Open the series page for your dataset on GEO (e.g. `GSE19804`) and note its platform ID (GSE19804 uses **GPL570**, the Affymetrix Human Genome U133 Plus 2.0 Array).
- Go to that platform's GEO page (`GPLxxxx`) and download its annotation table - GEO also mirrors these as plain files named like `GPL570-55999.txt`.
- Confirm the file has an ID column (probe IDs like `1007_s_at`) and a gene-symbol column (`Gene Symbol`) before use.

**2. The downloaded file often won't parse as-is.**
GEO platform annotation files are typically full tab-separated tables, but they:
- Start with several `#`-prefixed header/comment lines describing each column, *before* the real header row - these need to be stripped, or `pandas` will misread the whole file.
- Contain gene descriptions with embedded quotes, commas, and parentheses (e.g. multiple gene symbols separated by `///`) that can confuse a generic delimiter-autodetection parser and raise errors like `' ' expected after '"'`.

If you hit that error, open the file in a text editor first, count how many `#` lines precede the real header row, and strip them before passing it to `--id-map`, e.g. (PowerShell):

```powershell
Get-Content GPL570-55999.txt | Select-Object -Skip 16 | Set-Content GPL570_clean.txt
```

(replace `16` with however many comment lines your copy actually has), then run with `--id-map GPL570_clean.txt` instead of the raw download.

---

## Usage

### 1. Run the full pipeline

```bash
python deg_analysis.py --metadata meta.csv --expression expr.csv --outdir results
```

This loads the data, computes DEG statistics, saves the volcano plot and DEG report, runs the co-expression analysis, and - unless `--skip-go` is set - queries the Enrichr API for GO enrichment. The GO step can take a little while (a few seconds per gene list plus a short pause between calls), so budget for that on top of the rest, which typically finishes in seconds even on large matrices thanks to the vectorized statistics.

### 2. Customize column names, labels, and cutoffs

```bash
python deg_analysis.py --metadata meta.csv --expression expr.csv --sample-col sample_id --group-col condition --disease-label cancer --healthy-label control --fc-cutoff 1.0 --pval-cutoff 0.05 --corr-pval-cutoff 0.05 --outdir results
```

### 3. Run offline (skip GO enrichment)

```bash
python deg_analysis.py --metadata meta.csv --expression expr.csv --skip-go
```

### All available options

```bash
python deg_analysis.py --help
```

| Flag | Description | Default |
|---|---|---|
| `--metadata` | CSV file with sample metadata (sample IDs + group labels). **Required** | N/A |
| `--expression` | CSV file with the gene expression matrix (genes as rows, samples as columns). **Required** | N/A |
| `--sample-col` | Metadata column holding sample IDs | `sample_accession` |
| `--group-col` | Metadata column holding the condition label | `group` |
| `--disease-label` | Value in `--group-col` marking a diseased sample | `tumor` |
| `--healthy-label` | Value in `--group-col` marking a healthy/control sample | `normal` |
| `--id-map` | CSV/TSV mapping probe/feature IDs to gene symbols (e.g. a GEO platform annotation file) | none |
| `--id-map-id-col` | Column in `--id-map` holding the original probe/feature ID | `ID` |
| `--id-map-symbol-col` | Column in `--id-map` holding the gene symbol | `Gene Symbol` |
| `--id-map-agg` | How to aggregate multiple probes mapping to the same gene symbol (`mean`/`median`/`max`) | `mean` |
| `--fc-cutoff` | Absolute log2 fold-change cutoff for calling a gene a DEG | `1.0` |
| `--pval-cutoff` | Adjusted p-value cutoff for calling a gene a DEG | `0.05` |
| `--corr-pval-cutoff` | P-value cutoff for calling a co-expression pair significant | `0.05` |
| `--skip-go` | Skip the GO enrichment step (e.g. no internet access) | off |
| `--gene-sets` | gseapy/Enrichr gene-set library to use | `GO_Biological_Process_2021` |
| `--organism` | Organism passed to `gseapy.enrichr` | `human` |
| `--outdir` | Directory where all output files are written | `deg_output` |
| `-v`, `--verbose` | Enable verbose (DEBUG) logging | off |

---

## Statistical Methodology

**Fold change and significance.** For every gene, log2 fold change is computed from the group means (directly, if the data is already log-transformed; via a `log2` ratio with a small epsilon safeguard against division by zero, if the data is on a raw linear scale), and significance is assessed with a two-sample Welch's t-test (unequal variances assumed, no equal-variance assumption imposed on the two groups).

**Multiple-testing correction.** Raw p-values are adjusted with the Benjamini-Hochberg procedure, implemented directly with NumPy rather than pulled from a statistics package, to control the false discovery rate across all genes tested simultaneously.

**Classification.** A gene is called a DEG only if it clears both the fold-change cutoff (`--fc-cutoff`, log2 scale) and the adjusted p-value cutoff (`--pval-cutoff`) - genes that miss either threshold are labeled `equally-expressed`.

**Co-expression.** For every DEG pair (one up-regulated, one down-regulated), the Pearson correlation is computed across the diseased samples only, with the matching p-value obtained analytically from the correlation's t-distribution rather than from repeated calls to a hypothesis-testing function. Pairs are kept only if the correlation is statistically significant (`--corr-pval-cutoff`), and split into positive ("similar") and negative ("different") co-expression patterns.

**GO enrichment** is run as three independent Enrichr queries - up-regulated genes only, down-regulated genes only, and all DEGs together - so the biological processes implicated by each direction of change can be compared side by side.

---

## Project Structure

```
.
├── deg_analysis.py       # Main pipeline script (single entry point)
├── requirements.txt      # One-command dependency installation
├── README.md             # This file
└── deg_output/            # Generated: reports, plots, and GO enrichment results (created after a run)
```

---

## Output Files

Running the pipeline writes the following to `--outdir` (default `deg_output/`):

- `volcano_plot.png` - fold change vs. -log10 adjusted p-value, color-coded red (over-expressed) / blue (down-expressed) / grey (unchanged).
- `deg_report.txt` - tab-separated, one line per gene: `gene`, `fold_change`, `p_value`, `adj_p_value`, `t_stat`, `change_type`.
- `coexpression_similar.tsv` - significant positive-correlation up/down gene pairs, ranked descending by correlation.
- `coexpression_different.tsv` - significant negative-correlation up/down gene pairs, ranked ascending by correlation (most negative first).
- `go_enrichment_up/`, `go_enrichment_down/`, `go_enrichment_all/` - Enrichr GO Biological Process results for each gene subset (skipped entirely if `--skip-go` is set).

**Output example (console log):**

```
[INFO] Diseased matrix shape: (54675, 60) | Healthy matrix shape: (54675, 60)
[INFO] Detected that the expression matrix is in raw linear scale.
[INFO] Gene regulation breakdown:
change_type
equally-expressed    53812
over-expressed          421
down-expressed          442
[INFO] Volcano plot saved to deg_output/volcano_plot.png
[INFO] DEG report saved to deg_output/deg_report.txt
[INFO] Co-expression reports saved to deg_output/coexpression_similar.tsv and deg_output/coexpression_different.tsv
[INFO] Found 187 significant similar (positive) pairs and 203 significant different (negative) pairs.
[INFO] Done. All outputs written to /path/to/deg_output
```

---

## Troubleshooting

**`ValueError: --id-map file is missing required column(s)... Available columns: ['ID_REF', 'GSM...', ...]`**
You pointed `--id-map` at the expression matrix itself instead of a platform annotation file. `--id-map` needs a separate GPL annotation file with an ID column and a `Gene Symbol` column - see [Getting and cleaning a GPL annotation file](#getting-and-cleaning-a-gpl-annotation-file-common-pitfalls) above.

**`pandas.errors.ParserError: ' ' expected after '"'` while loading `--id-map`**
The raw GEO annotation file has leading `#` comment lines and/or quoted fields that trip up automatic delimiter detection. Strip the comment lines and re-save as a clean tab-separated file before passing it to `--id-map` (see the same section above for the exact command).

**GO enrichment logs "succeeded but returned 0 annotations" for every gene list**
Your expression matrix (and therefore `results["gene"]`) is still indexed by probe ID, not gene symbol - Enrichr can't match probe IDs to anything. Rerun with a working `--id-map` so genes are translated to symbols before the GO step.

---

## Notes and Limitations

- There is no default `--metadata`/`--expression` dataset bundled with this tool; both are always required.
- **GO enrichment requires gene symbols, not probe IDs.** If your expression matrix is indexed by microarray probe/feature ID (common for Affymetrix platforms), GO enrichment will silently return 0 annotations unless you map probes to gene symbols first with `--id-map` - see [A note on gene identifiers and GO enrichment](#a-note-on-gene-identifiers-and-go-enrichment). The pipeline detects likely probe IDs and warns if `--id-map` isn't supplied, and flags any Enrichr call that comes back with 0 annotations regardless.
- The script auto-detects whether the expression matrix is already log-transformed (max value < 25) and adjusts the fold-change formula accordingly - double-check this against your actual data if your values happen to fall near that boundary.
- The GO enrichment step depends on the public Enrichr API being reachable; it fails gracefully per gene list (logs a warning, moves on) rather than crashing the whole run, but you'll still need `--skip-go` for fully offline use.
- Co-expression significance uses an analytic p-value derived from the correlation's t-distribution rather than a permutation test - this is the standard approach for Pearson correlation but assumes approximately normally distributed residuals.
- This tool is intended for research and educational purposes and does **not** constitute clinical or diagnostic advice.

---
## Acknowledgements

This pipeline was originally developed as a term project for the **CIT672: Programming for Bioinformatics** course at **Nile University**. 

### Group (A) Members
1. **Akram Nader Salah**
2. **Mariam Mohamed Mahdy**
3. **Ala'a Ahmed Fathallah**
4. **Mohammed Fouad Mohammed**

**Project Topic:** *Differentially Expressed Genes (DEG) Analysis for Lung Cancer and Gene Ontology (GO) Enrichment*

`Thank you to everyone in the group for their collaboration, dedication, and shared efforts throughout this project.` 

A special thanks to our course instructor, **Dr. Ibrahim Mohamed Youssef**, for his invaluable guidance, effort, and support.

---
## License

This project is released under the [MIT License](https://opensource.org/licenses/MIT).
