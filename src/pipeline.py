"""
pipeline.py
===========
The scientific heart of the project. Each function is one logical stage of a
standard single-cell RNA-seq (scRNA-seq) workflow built on Scanpy:

    load_data            -> read the 10x count matrix into an AnnData object
    run_qc               -> compute QC metrics and remove low-quality cells
    normalize_and_select -> normalise, log-transform, pick variable genes, scale, PCA
    reduce_and_cluster   -> build the kNN graph, run UMAP and Leiden clustering
    find_markers         -> rank differentially expressed genes per cluster
    annotate_clusters    -> assign a cell type to each cluster from canonical markers

Every function takes and returns an AnnData object, so the stages compose like
a pipeline. The orchestration lives in `run_pipeline.py`.

AnnData primer
--------------
Scanpy stores everything in a single object, `adata`:
  * adata.X        -> the expression matrix (cells x genes)
  * adata.obs      -> per-cell metadata table (QC metrics, cluster labels, ...)
  * adata.var      -> per-gene metadata table (gene ids, "is highly variable", ...)
  * adata.obsm     -> per-cell matrices (PCA coordinates, UMAP coordinates, ...)
  * adata.uns      -> unstructured extras (results of rank_genes_groups, ...)
  * adata.raw      -> a frozen copy of the log-normalised data for all genes
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import scanpy as sc

from . import config as cfg


# --------------------------------------------------------------------------- #
# 0. Scanpy global settings
# --------------------------------------------------------------------------- #
def configure_scanpy() -> None:
    """Set verbosity, figure defaults and the global random seed."""
    sc.settings.verbosity = 1                 # 0=errors, 1=warnings, 2=info, 3=hints
    sc.settings.set_figure_params(dpi=150, dpi_save=200, facecolor="white")
    sc.settings.figdir = cfg.FIGURES_DIR      # where sc.pl.*(save=...) writes
    np.random.seed(cfg.RANDOM_SEED)


# --------------------------------------------------------------------------- #
# 1. Load
# --------------------------------------------------------------------------- #
def load_data() -> sc.AnnData:
    """
    Read the 10x Cell Ranger output (matrix.mtx + genes.tsv + barcodes.tsv)
    into an AnnData object.

    `var_names="gene_symbols"` uses human-readable gene names (CD3D, MS4A1, ...)
    rather than Ensembl IDs. `make_unique` appends suffixes if a symbol appears
    twice, which is required because some symbols are not unique.
    """
    adata = sc.read_10x_mtx(cfg.DATA_DIR, var_names="gene_symbols", cache=True)
    adata.var_names_make_unique()
    print(f"[load] {adata.n_obs} cells x {adata.n_vars} genes")
    return adata


# --------------------------------------------------------------------------- #
# 2. Quality control
# --------------------------------------------------------------------------- #
def run_qc(adata: sc.AnnData) -> sc.AnnData:
    """
    Compute QC metrics, visualise them, then filter out low-quality barcodes
    and rarely-detected genes.

    The three QC pillars in droplet scRNA-seq:
      * n_genes_by_counts : how many distinct genes a cell expresses.
                            Too few  -> empty droplet / debris.
                            Too many -> two cells captured together (doublet).
      * total_counts      : total transcripts (UMIs) per cell (sequencing depth).
      * pct_counts_mt     : fraction of counts from mitochondrial genes.
                            High  -> the cell membrane broke and cytoplasmic mRNA
                            leaked out, leaving mostly mitochondrial transcripts:
                            a hallmark of a stressed or dying cell.
    """
    # Basic gene/cell filtering first.
    sc.pp.filter_cells(adata, min_genes=cfg.MIN_GENES_PER_CELL)
    sc.pp.filter_genes(adata, min_cells=cfg.MIN_CELLS_PER_GENE)

    # Flag mitochondrial genes and let Scanpy compute the standard metrics.
    adata.var["mt"] = adata.var_names.str.startswith(cfg.MITO_PREFIX)
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
    )

    # --- figures: distributions BEFORE the cut --------------------------- #
    sc.pl.violin(
        adata,
        ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
        jitter=0.4, multi_panel=True, show=False, save="_qc_before.png",
    )
    sc.pl.scatter(
        adata, x="total_counts", y="pct_counts_mt",
        show=False, save="_qc_counts_vs_mito.png",
    )
    sc.pl.scatter(
        adata, x="total_counts", y="n_genes_by_counts",
        show=False, save="_qc_counts_vs_genes.png",
    )

    # --- apply the thresholds ------------------------------------------- #
    n_before = adata.n_obs
    adata = adata[adata.obs.n_genes_by_counts < cfg.MAX_GENES_PER_CELL, :].copy()
    adata = adata[adata.obs.pct_counts_mt < cfg.MAX_PCT_MITO, :].copy()
    print(f"[qc] kept {adata.n_obs}/{n_before} cells after filtering")

    sc.pl.violin(
        adata,
        ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
        jitter=0.4, multi_panel=True, show=False, save="_qc_after.png",
    )
    return adata


# --------------------------------------------------------------------------- #
# 3. Normalisation, feature selection, scaling, PCA
# --------------------------------------------------------------------------- #
def normalize_and_select(adata: sc.AnnData) -> sc.AnnData:
    """
    Turn raw integer counts into something comparable across cells, then reduce
    the gene space to the most informative dimensions.

    Steps:
      1. normalize_total : rescale each cell to TARGET_SUM counts (CP10K).
      2. log1p           : log(1 + x). Gene expression is highly skewed; the log
                           tames the dynamic range so a few very high counts do
                           not dominate everything downstream.
      3. highly_variable_genes : keep genes whose variance exceeds technical noise.
      4. adata.raw = adata     : freeze the full log-normalised matrix; we will
                                 plot markers from it later even after subsetting.
      5. subset to HVGs        : ~2k genes instead of ~14k -> faster, less noise.
      6. regress_out           : remove the linear effect of depth & mito fraction.
      7. scale                 : z-score each gene (mean 0, unit variance) so that
                                 PCA is not dominated by a handful of high-mean genes.
      8. pca                   : linear dimensionality reduction; the principal
                                 components become the input to the cell graph.
    """
    sc.pp.normalize_total(adata, target_sum=cfg.TARGET_SUM)
    sc.pp.log1p(adata)

    sc.pp.highly_variable_genes(
        adata,
        min_mean=cfg.HVG_MIN_MEAN,
        max_mean=cfg.HVG_MAX_MEAN,
        min_disp=cfg.HVG_MIN_DISP,
    )
    n_hvg = int(adata.var.highly_variable.sum())
    print(f"[hvg] selected {n_hvg} highly variable genes")
    sc.pl.highly_variable_genes(adata, show=False, save="_hvg.png")

    # Keep the complete log-normalised matrix for marker plotting / DE testing.
    adata.raw = adata

    # Work only with the informative genes from here on.
    adata = adata[:, adata.var.highly_variable].copy()

    sc.pp.regress_out(adata, cfg.REGRESS_OUT)
    sc.pp.scale(adata, max_value=cfg.SCALE_MAX_VALUE)

    sc.tl.pca(adata, svd_solver="arpack", random_state=cfg.RANDOM_SEED)
    sc.pl.pca_variance_ratio(
        adata, n_pcs=50, log=True, show=False, save="_pca_variance.png"
    )
    return adata


# --------------------------------------------------------------------------- #
# 4. Graph, UMAP, clustering
# --------------------------------------------------------------------------- #
def reduce_and_cluster(adata: sc.AnnData) -> sc.AnnData:
    """
    Build a k-nearest-neighbour graph in PCA space, embed it in 2D with UMAP for
    visualisation, and partition it into communities with the Leiden algorithm.

    Intuition: cells that are transcriptionally similar end up close together in
    PCA space. The kNN graph connects each cell to its k most similar cells.
    Leiden then finds densely connected groups in that graph - these groups are
    our candidate cell populations. UMAP is *only* for drawing the picture; the
    clustering is done on the high-dimensional graph, not on the 2D coordinates.
    """
    sc.pp.neighbors(
        adata, n_neighbors=cfg.N_NEIGHBORS, n_pcs=cfg.N_PCS,
        random_state=cfg.RANDOM_SEED,
    )
    sc.tl.umap(adata, random_state=cfg.RANDOM_SEED)
    sc.tl.leiden(
        adata,
        resolution=cfg.LEIDEN_RESOLUTION,
        random_state=cfg.RANDOM_SEED,
        flavor="igraph", n_iterations=2, directed=False,
    )
    n_clusters = adata.obs["leiden"].nunique()
    print(f"[cluster] Leiden found {n_clusters} clusters")

    sc.pl.umap(
        adata, color="leiden", legend_loc="on data", title="Leiden clusters",
        frameon=False, show=False, save="_clusters.png",
    )
    sc.pl.umap(
        adata, color=["n_genes_by_counts", "total_counts", "pct_counts_mt"],
        frameon=False, show=False, save="_qc_on_umap.png",
    )
    return adata


# --------------------------------------------------------------------------- #
# 5. Marker genes (differential expression per cluster)
# --------------------------------------------------------------------------- #
def find_markers(adata: sc.AnnData, n_top: int = 25) -> pd.DataFrame:
    """
    For every cluster, find the genes that are most specifically up-regulated
    relative to all other cells (a one-vs-rest Wilcoxon rank-sum test). These
    "marker genes" are what we use to decide what each cluster actually is.

    Returns a tidy DataFrame (cluster, rank, gene, logfoldchange, pval_adj, score)
    and also writes the top markers to results/marker_genes.csv.
    """
    sc.tl.rank_genes_groups(
        adata, groupby="leiden", method="wilcoxon",
        use_raw=True, n_genes=n_top,
    )
    sc.pl.rank_genes_groups(
        adata, n_genes=20, sharey=False, show=False, save="_marker_genes.png"
    )

    # Flatten Scanpy's structured result into a tidy table.
    res = adata.uns["rank_genes_groups"]
    groups = res["names"].dtype.names
    records = []
    for grp in groups:
        for rank in range(len(res["names"][grp])):
            records.append(
                {
                    "cluster": grp,
                    "rank": rank + 1,
                    "gene": res["names"][grp][rank],
                    "log2_fold_change": float(res["logfoldchanges"][grp][rank]),
                    "pval_adj": float(res["pvals_adj"][grp][rank]),
                    "score": float(res["scores"][grp][rank]),
                }
            )
    markers = pd.DataFrame.from_records(records)
    markers.to_csv(cfg.MARKER_CSV, index=False)
    print(f"[markers] wrote {cfg.MARKER_CSV.name} "
          f"({n_top} genes x {len(groups)} clusters)")
    return markers


# --------------------------------------------------------------------------- #
# 6. Automated, marker-based cell-type annotation
# --------------------------------------------------------------------------- #
def annotate_clusters(adata: sc.AnnData) -> sc.AnnData:
    """
    Assign a biological identity to each Leiden cluster.

    Rather than hard-coding "cluster 3 = CD8 T cells" (fragile, because cluster
    numbers shuffle between Scanpy versions), we score every cluster against the
    canonical marker sets in config.MARKER_SETS and take the best match:

      1. For each marker gene, take its mean log-normalised expression in each
         cluster (from adata.raw, which still has all genes).
      2. Z-score that gene across clusters, so a cluster is "high" for a gene
         relative to the other clusters, not in absolute terms.
      3. Average the z-scores within each marker set -> one score per (cluster,
         cell type).
      4. Assign each cluster the cell type with the highest score (argmax).

    A dotplot and the final annotated UMAP are saved for visual confirmation.
    """
    # Pull the full log-normalised matrix back out of adata.raw as a DataFrame.
    raw = adata.raw.to_adata()
    available = {g for g in raw.var_names}
    clusters = sorted(adata.obs["leiden"].cat.categories, key=int)

    # mean expression of every marker gene in every cluster
    all_markers = sorted({g for genes in cfg.MARKER_SETS.values()
                          for g in genes if g in available})
    expr = pd.DataFrame(index=clusters, columns=all_markers, dtype=float)
    for cl in clusters:
        mask = (adata.obs["leiden"] == cl).values
        sub = raw[mask, all_markers].X
        sub = np.asarray(sub.mean(axis=0)).ravel()
        expr.loc[cl] = sub

    # z-score each gene across clusters
    zexpr = (expr - expr.mean(axis=0)) / (expr.std(axis=0) + 1e-9)

    # score each cluster against each cell-type marker set
    score = pd.DataFrame(index=clusters, columns=list(cfg.MARKER_SETS), dtype=float)
    for ctype, genes in cfg.MARKER_SETS.items():
        present = [g for g in genes if g in available]
        score[ctype] = zexpr[present].mean(axis=1)

    mapping = score.idxmax(axis=1).to_dict()           # cluster -> cell type
    best_score = score.max(axis=1).to_dict()

    # Save the decision table so the annotation is transparent / auditable.
    mapping_df = (
        score.assign(assigned=pd.Series(mapping), confidence=pd.Series(best_score))
        .rename_axis("cluster").reset_index()
    )
    mapping_df.to_csv(cfg.MAPPING_CSV, index=False)
    print("[annotate] cluster -> cell type:")
    for cl in clusters:
        print(f"           cluster {cl:>2} -> {mapping[cl]:<18} "
              f"(score {best_score[cl]:+.2f})")

    # Attach the labels to every cell.
    adata.obs["cell_type"] = (
        adata.obs["leiden"].map(mapping).astype("category")
    )

    # --- figures --------------------------------------------------------- #
    # dotplot of canonical markers grouped by cluster (visual sanity check)
    dot_genes = {ct: [g for g in genes if g in available]
                 for ct, genes in cfg.MARKER_SETS.items()}
    sc.pl.dotplot(
        adata, dot_genes, groupby="leiden", use_raw=True,
        standard_scale="var", show=False, save="_marker_dotplot.png",
    )
    # the headline figure: UMAP coloured by annotated cell type
    sc.pl.umap(
        adata, color="cell_type", legend_loc="right margin",
        title="PBMC 3k - annotated cell types", frameon=False,
        show=False, save="_celltypes.png",
    )
    # feature plots: one canonical marker per lineage on the UMAP
    spotlight = [g for g in cfg.SPOTLIGHT_MARKERS if g in available]
    sc.pl.umap(
        adata, color=spotlight, use_raw=True, frameon=False,
        ncols=3, show=False, save="_marker_featureplots.png",
    )

    # --- cell-type proportions ------------------------------------------ #
    counts = adata.obs["cell_type"].value_counts()
    counts.rename_axis("cell_type").reset_index(name="n_cells").to_csv(
        cfg.COUNTS_CSV, index=False
    )

    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7, 4))
    counts.sort_values().plot(kind="barh", ax=ax, color="#4C72B0")
    ax.set_xlabel("number of cells")
    ax.set_ylabel("")
    ax.set_title("Cell-type composition of PBMC 3k")
    for i, v in enumerate(counts.sort_values().values):
        ax.text(v + 5, i, str(int(v)), va="center", fontsize=9)
    fig.tight_layout()
    return adata


# --------------------------------------------------------------------------- #
# 7. Housekeeping: give the saved figures clean, ordered names
# --------------------------------------------------------------------------- #
def tidy_figures() -> None:
    """Rename Scanpy's auto-named PNGs to the clean sequence in config.FIGURE_RENAMES."""
    for src_name, dst_name in cfg.FIGURE_RENAMES.items():
        src = cfg.FIGURES_DIR / src_name
        if src.exists():
            src.replace(cfg.FIGURES_DIR / dst_name)
    print(f"[tidy] renamed figures -> {len(cfg.FIGURE_RENAMES)} clean filenames")
