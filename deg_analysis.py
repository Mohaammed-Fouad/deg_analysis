#!/usr/bin/env python3
"""
deg_analysis.py

Command-line tool for Differentially Expressed Gene (DEG) analysis.

Given a gene-expression matrix and a sample metadata file that labels each
sample as "diseased" or "healthy", this script:

  1. Loads and splits the expression matrix into disease / healthy groups.
     Optionally maps probe/feature IDs (e.g. Affymetrix probe IDs like
     "1007_s_at") to gene symbols via --id-map, since GO enrichment
     (Enrichr) only matches on gene symbols -- see "A note on gene
     identifiers and GO enrichment" below.
  2. Computes log2 fold change, Welch's t-test statistic, and p-value for
     every gene (no DEG-calling packages are used -- only NumPy/SciPy).
  3. Adjusts p-values for multiple testing with the Benjamini-Hochberg
     procedure (implemented from scratch).
  4. Classifies every gene as over-expressed, down-expressed, or
     equally-expressed.
  5. Draws a volcano plot (matplotlib only).
  6. Saves a tab-separated DEG report (gene, fold_change, p_value,
     adj_p_value, t_stat, change_type).
  7. Runs a pairwise co-expression analysis between every up-regulated gene
     and every down-regulated gene (Pearson correlation on the diseased
     samples), and reports statistically-significant similar (positive) and
     different (negative) co-expression pairs, ranked by correlation.
  8. Runs GO Biological Process enrichment (via gseapy/Enrichr) separately
     for the up-regulated genes, the down-regulated genes, and all DEGs
     together.

This is a straight command-line port of the original Jupyter-notebook
prototype: the statistics, the BH correction, the classification rules, the
volcano-plot mechanics, and the co-expression logic are all unchanged. What
changed is packaging: everything now lives in importable functions driven
by a proper CLI (argparse) with a main() entry point, so the whole pipeline
can be run non-interactively, with logging and configurable parameters,
instead of by hand cell-by-cell.

A note on gene identifiers and GO enrichment
---------------------------------------------
Enrichr (used for the GO enrichment step) matches submitted genes by
*gene symbol* (e.g. "TP53", "EGFR"). Many microarray expression matrices
are indexed by *probe/feature ID* instead (e.g. Affymetrix IDs like
"1007_s_at"), which will never match anything in Enrichr's gene sets. This
fails silently: the API call still "succeeds," it just returns an empty
results table (0 annotations) for every gene list. If your gene index looks
like probe IDs, use --id-map to supply a CSV that maps each probe/feature
ID to its gene symbol (e.g. a GEO platform annotation file, which for many
Affymetrix platforms already uses "ID" and "Gene Symbol" as column names --
the defaults for --id-map-id-col / --id-map-symbol-col). The script also
logs a warning if it detects Affymetrix-style probe IDs in --id-map's
absence.

Example
-------
    python deg_analysis.py --metadata meta.csv --expression expr.csv --outdir results

    python deg_analysis.py --metadata meta.csv --expression expr.csv --disease-label tumor --healthy-label normal --fc-cutoff 1.0 --pval-cutoff 0.05 --skip-go

    python deg_analysis.py --metadata meta.csv --expression expr.csv --id-map GPL570-55999.txt --id-map-id-col ID --id-map-symbol-col "Gene Symbol"
"""


"""
Acknowledgements
----------------
This pipeline was originally developed as a term project for the CIT672: 
Programming for Bioinformatics course at Nile University. 

Group (A) members:
    1. Akram Nader Salah
    2. Mariam Mohamed Mahdy
    3. Ala'a Ahmed Fathallah
    4. Mohammed Fouad Mohammed

Topic: 
    Differentially Expressed Genes (DEG) Analysis for Lung Cancer and Gene Ontology (GO) Enrichment

Thank you to everyone in the group for their collaboration and effort throughout this project. 

A special thanks to our instructor, Dr. Ibrahim Mohamed Youssef, for his effort and support during this course.
"""


import argparse
import logging
import os
import re
import sys
import time

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")  # safe for headless / CLI use
import matplotlib.pyplot as plt

logger = logging.getLogger("deg_analysis")

# Affymetrix-style probe IDs look like "1007_s_at", "121_at", "AFFX-BioB-5_at".
# If most of the gene index matches this pattern and no --id-map was given,
# we warn the user, since Enrichr will not be able to match these to any
# GO gene set.
_PROBE_ID_PATTERN = re.compile(r"^(AFFX-)?[\w.-]+_(a_|s_|x_|i_)?at$", re.IGNORECASE)




# --------------------------------------------------------------------------- #
# Step 1: Load and split the data
# --------------------------------------------------------------------------- #
# Gene identifier mapping (probe/feature ID -> gene symbol)
# --------------------------------------------------------------------------- #
def looks_like_probe_ids(gene_index, sample_size=100):
    """Heuristically detect Affymetrix-style probe IDs (e.g. '1007_s_at').

    Enrichr matches genes by symbol, not by probe ID, so if most of the
    index looks like this, GO enrichment will silently return 0
    annotations unless the caller supplies --id-map.
    """
    sample = [str(g) for g in list(gene_index)[:sample_size]]
    if not sample:
        return False
    matches = sum(1 for g in sample if _PROBE_ID_PATTERN.match(g))
    return (matches / len(sample)) > 0.5


def load_id_mapping(id_map_path, id_col, symbol_col):
    """Load a probe/feature-ID -> gene-symbol CSV (e.g. a GEO platform
    annotation file) into a {id: symbol} dict."""
    logger.info("Loading gene ID mapping from %s", id_map_path)
    map_df = pd.read_csv(id_map_path, sep=None, engine="python")

    if id_col not in map_df.columns or symbol_col not in map_df.columns:
        raise ValueError(
            f"--id-map file is missing required column(s). Looked for "
            f"id column '{id_col}' and symbol column '{symbol_col}'. "
            f"Available columns: {list(map_df.columns)}"
        )

    map_df = map_df[[id_col, symbol_col]].dropna()
    map_df[id_col] = map_df[id_col].astype(str)
    map_df[symbol_col] = map_df[symbol_col].astype(str).str.strip()
    map_df = map_df[map_df[symbol_col] != ""]

    mapping = dict(zip(map_df[id_col], map_df[symbol_col]))
    logger.info("Loaded %d probe/feature ID -> gene symbol mappings.", len(mapping))
    return mapping


def map_to_gene_symbols(expression_df, mapping, agg="mean"):
    """Translate an expression matrix's row index from probe/feature IDs to
    gene symbols using `mapping`, dropping unmapped rows and collapsing
    multiple probes that map to the same symbol via `agg`."""
    original_n = expression_df.shape[0]

    mapped_index = expression_df.index.to_series().astype(str).map(mapping)
    keep_mask = mapped_index.notna() & (mapped_index.str.strip() != "")

    mapped_df = expression_df.loc[keep_mask.values].copy()
    mapped_df.index = mapped_index[keep_mask].values

    dropped = original_n - mapped_df.shape[0]
    if dropped:
        logger.warning("Dropped %d/%d rows with no gene symbol in --id-map.",
                        dropped, original_n)

    n_before_collapse = mapped_df.shape[0]
    collapsed_df = mapped_df.groupby(mapped_df.index).agg(agg)
    n_after_collapse = collapsed_df.shape[0]
    if n_after_collapse < n_before_collapse:
        logger.info(
            "Collapsed %d probes mapping to duplicate gene symbols into %d "
            "unique genes using '%s' aggregation.",
            n_before_collapse, n_after_collapse, agg,
        )

    if collapsed_df.shape[0] == 0:
        raise ValueError(
            "No rows survived gene-symbol mapping. Double-check --id-map-id-col "
            "and --id-map-symbol-col against your --id-map file's actual columns, "
            "and that the IDs it contains actually match the expression matrix's "
            "gene index."
        )

    return collapsed_df


# --------------------------------------------------------------------------- #
# Step 1: Load and split the data
# --------------------------------------------------------------------------- #
def load_data(metadata_path, expression_path, sample_col, group_col,
              disease_label, healthy_label, id_map_path=None,
              id_map_id_col="ID", id_map_symbol_col="Gene Symbol",
              id_map_agg="mean"):
    """Load the metadata + expression matrix and split into disease/healthy.

    If `id_map_path` is given, the expression matrix's gene index is first
    translated from probe/feature IDs to gene symbols (see the module
    docstring's note on gene identifiers and GO enrichment). Otherwise, the
    gene index is left as-is, but a warning is logged if it looks like
    Affymetrix-style probe IDs, since those will not match anything in
    Enrichr's GO gene sets later on.
    """
    logger.info("Loading metadata from %s", metadata_path)
    metadata_df = pd.read_csv(metadata_path)

    if group_col not in metadata_df.columns:
        raise ValueError(
            f"Group column '{group_col}' not found in metadata. "
            f"Available columns: {list(metadata_df.columns)}"
        )
    if sample_col not in metadata_df.columns:
        raise ValueError(
            f"Sample-ID column '{sample_col}' not found in metadata. "
            f"Available columns: {list(metadata_df.columns)}"
        )

    disease_sample_ids = metadata_df[metadata_df[group_col] == disease_label][sample_col].tolist()
    healthy_sample_ids = metadata_df[metadata_df[group_col] == healthy_label][sample_col].tolist()

    logger.info("Found %d disease samples and %d healthy samples.",
                len(disease_sample_ids), len(healthy_sample_ids))

    if len(disease_sample_ids) == 0 or len(healthy_sample_ids) == 0:
        raise ValueError(
            "One of the groups is empty. Check --disease-label / --healthy-label "
            f"against the values actually present in column '{group_col}'."
        )

    logger.info("Loading expression matrix from %s", expression_path)
    expression_df = pd.read_csv(expression_path, index_col=0)

    if id_map_path:
        mapping = load_id_mapping(id_map_path, id_map_id_col, id_map_symbol_col)
        expression_df = map_to_gene_symbols(expression_df, mapping, agg=id_map_agg)
    elif looks_like_probe_ids(expression_df.index):
        logger.warning(
            "The gene index looks like Affymetrix-style probe IDs (e.g. "
            "'1007_s_at') rather than gene symbols. GO enrichment (Enrichr) "
            "matches by gene symbol, so it will silently return 0 "
            "annotations for every list unless you map probes to symbols "
            "first with --id-map (e.g. a GEO platform annotation file)."
        )

    missing_disease = [s for s in disease_sample_ids if s not in expression_df.columns]
    missing_healthy = [s for s in healthy_sample_ids if s not in expression_df.columns]
    if missing_disease or missing_healthy:
        raise ValueError(
            "Some sample IDs listed in the metadata are missing from the "
            f"expression matrix columns. Missing disease samples: {missing_disease[:5]}"
            f"{'...' if len(missing_disease) > 5 else ''}; "
            f"missing healthy samples: {missing_healthy[:5]}"
            f"{'...' if len(missing_healthy) > 5 else ''}"
        )

    disease_df = expression_df[disease_sample_ids]
    healthy_df = expression_df[healthy_sample_ids]

    logger.info("Diseased matrix shape: %s | Healthy matrix shape: %s",
                disease_df.shape, healthy_df.shape)

    return expression_df, disease_df, healthy_df


# --------------------------------------------------------------------------- #
# Step 2: Fold change + t-test per gene
# --------------------------------------------------------------------------- #
def calculate_stats(diseased_df, healthy_df, epsilon=1e-5):
    """Compute log2 fold change, t-statistic, and p-value for every gene.

    Fully vectorized: means, fold changes, and Welch's t-test are all
    computed across every gene at once (axis=1 array operations) instead of
    looping gene-by-gene. This is functionally equivalent to the original
    per-gene loop, but scales far better on large genomic datasets.
    """
    gene_names = diseased_df.index.to_numpy()
    diseased_vals = diseased_df.to_numpy(dtype=float)  # shape: (n_genes, n_disease_samples)
    healthy_vals = healthy_df.to_numpy(dtype=float)    # shape: (n_genes, n_healthy_samples)

    # Detect if data is already log-transformed by checking maximum values
    # (if max value is small, e.g., < 25, it's already log-transformed)
    is_already_log = diseased_vals.max() < 25

    if is_already_log:
        logger.info("Detected that the expression matrix is ALREADY log-transformed.")
    else:
        logger.info("Detected that the expression matrix is in raw linear scale.")

    diseased_mean = diseased_vals.mean(axis=1)
    healthy_mean = healthy_vals.mean(axis=1)

    # Calculate log2 fold change based on data scale, across all genes at once
    if is_already_log:
        # If already logged, log(A) - log(B) = log(A/B)
        log2_fc = diseased_mean - healthy_mean
    else:
        # Raw linear numbers: add a small epsilon to both means so a zero
        # mean in either group can never trigger a division-by-zero or
        # log2(0) error.
        log2_fc = np.log2((diseased_mean + epsilon) / (healthy_mean + epsilon))

    # Run Welch's t-test for every gene simultaneously (axis=1 = across samples)
    t_stat, p_value = stats.ttest_ind(diseased_vals, healthy_vals, axis=1, equal_var=False)

    # Cleanly mask any NaN results (e.g. from zero-variance rows)
    nan_mask = np.isnan(p_value)
    p_value = np.where(nan_mask, 1.0, p_value)
    t_stat = np.where(nan_mask, 0.0, t_stat)

    return pd.DataFrame({
        "gene": gene_names,
        "fold_change": log2_fc,
        "p_value": p_value,
        "t_stat": t_stat,
    })


# --------------------------------------------------------------------------- #
# Step 3: Benjamini-Hochberg p-value adjustment (implemented from scratch)
# --------------------------------------------------------------------------- #
def adjust_pvalues(p_values):
    """Benjamini-Hochberg FDR correction, fully vectorized with NumPy
    (no external multiple-testing package used).

    Same algorithm as before: sort ascending, scale each p-value by
    n / rank, enforce monotonicity walking from the largest rank back to
    the smallest, clip to [0, 1], and restore the original gene order —
    just expressed as array operations instead of explicit Python loops.
    """
    p_values = np.asarray(p_values, dtype=float)
    n = p_values.shape[0]

    # sort by p-value, smallest first, and remember the original positions
    order = np.argsort(p_values)
    ranked = p_values[order]
    ranks = np.arange(1, n + 1)

    # scale each sorted p-value by n / rank
    adjusted_sorted = ranked * n / ranks

    # enforce the "never decreasing going backwards" rule via a reversed
    # running minimum (equivalent to the original backward loop)
    adjusted_sorted = np.minimum.accumulate(adjusted_sorted[::-1])[::-1]

    # cap values at 1.0 (a p-value can't be more than 1)
    adjusted_sorted = np.clip(adjusted_sorted, 0.0, 1.0)

    # put everything back in the original gene order
    final_adjusted = np.empty(n, dtype=float)
    final_adjusted[order] = adjusted_sorted

    return final_adjusted.tolist()


# --------------------------------------------------------------------------- #
# Step 4: Classify genes
# --------------------------------------------------------------------------- #
def classify_genes(results, fold_change_cutoff=1.0, p_value_cutoff=0.05):
    """Label each gene as over-expressed, down-expressed, or equally-expressed.

    Fully vectorized with np.select: builds one boolean condition per label
    across the whole DataFrame at once, instead of looping row by row.
    Function signature and return type are unchanged.
    """
    fold_change = results["fold_change"].to_numpy()
    adj_p = results["adj_p_value"].to_numpy()

    conditions = [
        (adj_p < p_value_cutoff) & (fold_change >= fold_change_cutoff),
        (adj_p < p_value_cutoff) & (fold_change <= -fold_change_cutoff),
    ]
    choices = ["over-expressed", "down-expressed"]

    results["change_type"] = np.select(conditions, choices, default="equally-expressed")
    return results


# --------------------------------------------------------------------------- #
# Step 5: Volcano plot
# --------------------------------------------------------------------------- #
def make_volcano_plot(results, outpath, fold_change_cutoff=1.0, p_value_cutoff=0.05):
    """Draw and save the volcano plot (fold change vs -log10 adjusted p)."""
    colors = []
    for change_type in results["change_type"]:
        if change_type == "over-expressed":
            colors.append("red")
        elif change_type == "down-expressed":
            colors.append("blue")
        else:
            colors.append("grey")

    y_values = -np.log10(results["adj_p_value"].replace(0, np.nextafter(0, 1)))

    plt.figure(figsize=(8, 6))
    plt.scatter(results["fold_change"], y_values, c=colors, alpha=0.6)

    plt.axhline(-np.log10(p_value_cutoff), color="black", linestyle="--")
    plt.axvline(fold_change_cutoff, color="black", linestyle="--")
    plt.axvline(-fold_change_cutoff, color="black", linestyle="--")

    plt.xlabel("log2 Fold Change")
    plt.ylabel("-log10 Adjusted p-value")
    plt.title("Volcano Plot: Diseased vs Healthy")
    plt.savefig(outpath)
    plt.close()
    logger.info("Volcano plot saved to %s", outpath)


# --------------------------------------------------------------------------- #
# Step 6: Save the DEG report
# --------------------------------------------------------------------------- #
def save_report(results, filename):
    """Write one line per gene: name, fold change, p-value, adj p-value,
    t-stat, change type (tab-separated).

    Uses pandas' CSV writer instead of manual string concatenation. Column
    order and header text are unchanged, fold_change and t_stat are rounded
    to 4 decimal places (exactly as before), and no extra index column is
    written.
    """
    export_columns = ["gene", "fold_change", "p_value", "adj_p_value", "t_stat", "change_type"]
    export_df = results[export_columns].copy()
    export_df["fold_change"] = export_df["fold_change"].round(4)
    export_df["t_stat"] = export_df["t_stat"].round(4)

    export_df.to_csv(filename, sep="\t", index=False, header=True)

    logger.info("DEG report saved to %s", filename)


# --------------------------------------------------------------------------- #
# Step 7: Co-expression analysis (up genes vs down genes)
# --------------------------------------------------------------------------- #
def coexpression_analysis(diseased_df, results, corr_pvalue_cutoff=0.05):
    """Pearson correlation between every up-regulated and down-regulated gene
    (computed on the diseased samples). Returns two DataFrames: significant
    positive ("similar") and significant negative ("different") pairs, each
    ranked by correlation strength.

    Fully vectorized: instead of looping over every (up_gene, down_gene)
    pair and calling `scipy.stats.pearsonr` one pair at a time, this builds
    the whole up-genes x down-genes correlation matrix in one shot with
    `np.corrcoef`, then derives the matching p-values analytically via the
    standard t-distribution transformation of the Pearson correlation
    coefficient:

        t = r * sqrt((n - 2) / (1 - r**2))
        p = scipy.stats.t.sf(|t|, df=n - 2) * 2

    where n is the number of diseased samples. This is mathematically
    equivalent to calling `pearsonr` for every pair, but computes all pairs
    at once.
    """
    up_genes = results[results["change_type"] == "over-expressed"]["gene"].tolist()
    down_genes = results[results["change_type"] == "down-expressed"]["gene"].tolist()

    empty_cols = ["up_gene", "down_gene", "correlation", "p_value"]
    if len(up_genes) == 0 or len(down_genes) == 0:
        return pd.DataFrame(columns=empty_cols), pd.DataFrame(columns=empty_cols)

    up_matrix = diseased_df.loc[up_genes].to_numpy(dtype=float)      # (n_up, n_samples)
    down_matrix = diseased_df.loc[down_genes].to_numpy(dtype=float)  # (n_down, n_samples)
    n_up, n_down = up_matrix.shape[0], down_matrix.shape[0]
    n_samples = diseased_df.shape[1]
    df = n_samples - 2

    # Correlate every up gene against every down gene in one call: stack
    # both matrices, take the full correlation matrix, then slice out the
    # up-vs-down block.
    full_corr = np.corrcoef(np.vstack([up_matrix, down_matrix]))
    corr_matrix = full_corr[:n_up, n_up:]  # shape: (n_up, n_down)

    with np.errstate(divide="ignore", invalid="ignore"):
        t_stat = corr_matrix * np.sqrt(df / (1 - corr_matrix ** 2))
        p_matrix = stats.t.sf(np.abs(t_stat), df=df) * 2

    # Flatten into one row per (up_gene, down_gene) pair. NaN correlations
    # (e.g. from a zero-variance gene) produce NaN p-values, which safely
    # fail the p < cutoff comparison below and are dropped automatically.
    up_idx, down_idx = np.meshgrid(np.arange(n_up), np.arange(n_down), indexing="ij")
    flat_up = np.array(up_genes)[up_idx.ravel()]
    flat_down = np.array(down_genes)[down_idx.ravel()]
    flat_corr = corr_matrix.ravel()
    flat_p = p_matrix.ravel()

    significant_mask = flat_p < corr_pvalue_cutoff

    positive_mask = significant_mask & (flat_corr > 0)
    negative_mask = significant_mask & (flat_corr < 0)

    similar_df = pd.DataFrame({
        "up_gene": flat_up[positive_mask],
        "down_gene": flat_down[positive_mask],
        "correlation": flat_corr[positive_mask],
        "p_value": flat_p[positive_mask],
    })
    different_df = pd.DataFrame({
        "up_gene": flat_up[negative_mask],
        "down_gene": flat_down[negative_mask],
        "correlation": flat_corr[negative_mask],
        "p_value": flat_p[negative_mask],
    })

    if len(similar_df) > 0:
        similar_df = similar_df.sort_values("correlation", ascending=False)
    if len(different_df) > 0:
        different_df = different_df.sort_values("correlation", ascending=True)

    return similar_df, different_df


def save_coexpression_reports(similar_df, different_df, outdir):
    """Save the similar / different co-expression pair tables to disk."""
    similar_path = os.path.join(outdir, "coexpression_similar.tsv")
    different_path = os.path.join(outdir, "coexpression_different.tsv")

    if len(similar_df) > 0:
        similar_df.to_csv(similar_path, sep="\t", index=False)
    else:
        pd.DataFrame(columns=["up_gene", "down_gene", "correlation", "p_value"]).to_csv(
            similar_path, sep="\t", index=False)

    if len(different_df) > 0:
        different_df.to_csv(different_path, sep="\t", index=False)
    else:
        pd.DataFrame(columns=["up_gene", "down_gene", "correlation", "p_value"]).to_csv(
            different_path, sep="\t", index=False)

    logger.info("Co-expression reports saved to %s and %s", similar_path, different_path)


# --------------------------------------------------------------------------- #
# Step 8: GO enrichment
# --------------------------------------------------------------------------- #
def run_go_enrichment(results, outdir, gene_sets, organism, sleep_seconds=3):
    """Run GO Biological Process enrichment for up-only, down-only, and all DEGs."""
    import gseapy as gp  # imported lazily so --skip-go doesn't require it installed

    up_genes = results[results["change_type"] == "over-expressed"]["gene"].tolist()
    down_genes = results[results["change_type"] == "down-expressed"]["gene"].tolist()
    all_degs = up_genes + down_genes

    logger.info("Sending profiles to Enrichr API: %d up, %d down genes...",
                len(up_genes), len(down_genes))

    for label, gene_list in [("up", up_genes), ("down", down_genes), ("all", all_degs)]:
        if len(gene_list) == 0:
            logger.info("%s -> Skipped: gene list is empty.", label)
            continue
        try:
            enr = gp.enrichr(
                gene_list=gene_list,
                gene_sets=gene_sets,
                organism=organism,
                outdir=os.path.join(outdir, f"go_enrichment_{label}"),
                no_plot=True,
            )
            n_annotations = len(enr.results)
            if n_annotations == 0:
                logger.warning(
                    "%s -> Enrichr call succeeded but returned 0 annotations. "
                    "This almost always means the submitted gene identifiers "
                    "don't match Enrichr's gene symbols (e.g. probe IDs like "
                    "'1007_s_at' instead of symbols like 'TP53'). If your "
                    "expression matrix is indexed by probe/feature ID, rerun "
                    "with --id-map to translate to gene symbols first.",
                    label,
                )
            else:
                logger.info("Success: %s pathway processing saved. Found %d annotations.",
                            label, n_annotations)
        except Exception as e:
            logger.warning("Error processing %s: %s", label, e)
        time.sleep(sleep_seconds)


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def build_arg_parser():
    parser = argparse.ArgumentParser(
        description="Identify differentially expressed genes (DEGs) between a "
                    "diseased group and a healthy group, and run downstream "
                    "co-expression and GO enrichment analyses.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    io_group = parser.add_argument_group("Input files")
    io_group.add_argument("--metadata", required=True,
                           help="CSV file with sample metadata (sample IDs + group labels).")
    io_group.add_argument("--expression", required=True,
                           help="CSV file with the gene expression matrix "
                                "(genes as rows, samples as columns; first column = gene name).")
    io_group.add_argument("--sample-col", default="sample_accession",
                           help="Column in the metadata file holding sample IDs.")
    io_group.add_argument("--group-col", default="group",
                           help="Column in the metadata file holding the condition label.")
    io_group.add_argument("--disease-label", default="tumor",
                           help="Value in --group-col that marks a diseased sample.")
    io_group.add_argument("--healthy-label", default="normal",
                           help="Value in --group-col that marks a healthy/control sample.")
    io_group.add_argument("--id-map", default=None,
                           help="Optional CSV mapping probe/feature IDs to gene symbols "
                                "(e.g. a GEO platform annotation file). Recommended for "
                                "microarray platforms such as Affymetrix, whose probe IDs "
                                "(e.g. '1007_s_at') won't match anything in Enrichr's GO "
                                "gene sets otherwise -- see the module docstring.")
    io_group.add_argument("--id-map-id-col", default="ID",
                           help="Column in --id-map holding the original probe/feature ID.")
    io_group.add_argument("--id-map-symbol-col", default="Gene Symbol",
                           help="Column in --id-map holding the gene symbol to map to.")
    io_group.add_argument("--id-map-agg", default="mean", choices=["mean", "median", "max"],
                           help="How to aggregate multiple probes that map to the same "
                                "gene symbol.")

    stats_group = parser.add_argument_group("DEG calling thresholds")
    stats_group.add_argument("--fc-cutoff", type=float, default=1.0,
                              help="Absolute log2 fold-change cutoff for calling a gene a DEG.")
    stats_group.add_argument("--pval-cutoff", type=float, default=0.05,
                              help="Adjusted p-value cutoff for calling a gene a DEG.")
    stats_group.add_argument("--corr-pval-cutoff", type=float, default=0.05,
                              help="P-value cutoff for calling a co-expression pair significant.")

    go_group = parser.add_argument_group("GO enrichment")
    go_group.add_argument("--skip-go", action="store_true",
                           help="Skip the GO enrichment step (e.g. no internet access).")
    go_group.add_argument("--gene-sets", default="GO_Biological_Process_2021",
                           help="gseapy/Enrichr gene-set library to use.")
    go_group.add_argument("--organism", default="human",
                           help="Organism passed to gseapy.enrichr.")

    out_group = parser.add_argument_group("Output")
    out_group.add_argument("--outdir", default="deg_output",
                            help="Directory where all output files are written.")

    parser.add_argument("-v", "--verbose", action="store_true",
                         help="Enable verbose (DEBUG) logging.")

    return parser


def main(argv=None):
    parser = build_arg_parser()
    args = parser.parse_args(argv)

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%H:%M:%S",
    )
    # matplotlib's font-matching logs are extremely noisy at DEBUG and add
    # no value here, so keep them quiet even when --verbose is set.
    logging.getLogger("matplotlib").setLevel(logging.WARNING)

    os.makedirs(args.outdir, exist_ok=True)

    # Step 1: load data
    _, disease_df, healthy_df = load_data(
        args.metadata, args.expression, args.sample_col, args.group_col,
        args.disease_label, args.healthy_label,
        id_map_path=args.id_map, id_map_id_col=args.id_map_id_col,
        id_map_symbol_col=args.id_map_symbol_col, id_map_agg=args.id_map_agg,
    )

    # Step 2: fold change + t-test
    results = calculate_stats(disease_df, healthy_df)

    # Step 3: adjust p-values
    results["adj_p_value"] = adjust_pvalues(results["p_value"].tolist())

    # Step 4: classify genes
    results = classify_genes(results, fold_change_cutoff=args.fc_cutoff,
                              p_value_cutoff=args.pval_cutoff)

    counts = results["change_type"].value_counts()
    logger.info("Gene regulation breakdown:\n%s", counts.to_string())

    # Step 5: volcano plot
    volcano_path = os.path.join(args.outdir, "volcano_plot.png")
    make_volcano_plot(results, volcano_path, fold_change_cutoff=args.fc_cutoff,
                       p_value_cutoff=args.pval_cutoff)

    # Step 6: DEG report
    report_path = os.path.join(args.outdir, "deg_report.txt")
    save_report(results, filename=report_path)

    # Step 7: co-expression analysis
    similar_df, different_df = coexpression_analysis(
        disease_df, results, corr_pvalue_cutoff=args.corr_pval_cutoff)
    save_coexpression_reports(similar_df, different_df, args.outdir)
    logger.info("Found %d significant similar (positive) pairs and %d "
                "significant different (negative) pairs.",
                len(similar_df), len(different_df))

    # Step 8: GO enrichment
    if args.skip_go:
        logger.info("Skipping GO enrichment (--skip-go).")
    else:
        run_go_enrichment(results, args.outdir, gene_sets=[args.gene_sets],
                           organism=args.organism)

    logger.info("Done. All outputs written to %s", os.path.abspath(args.outdir))


if __name__ == "__main__":
    sys.exit(main())
