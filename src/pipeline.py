"""Analysis steps for the PBMC 3k dataset, one function per stage.

Each function takes and returns an AnnData object; run_pipeline.py calls them
in order.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import scanpy as sc

from . import config as cfg


def configure_scanpy() -> None:
    """Verbosity, figure defaults and the global seed."""
    sc.settings.verbosity = 1
    sc.settings.set_figure_params(dpi=150, dpi_save=200, facecolor="white")
    sc.settings.figdir = cfg.FIGURES_DIR
    np.random.seed(cfg.RANDOM_SEED)


def load_data() -> sc.AnnData:
    """Read the 10x matrix into AnnData, using gene symbols as var_names."""
    adata = sc.read_10x_mtx(cfg.DATA_DIR, var_names="gene_symbols", cache=True)
    adata.var_names_make_unique()
    print(f"[load] {adata.n_obs} cells x {adata.n_vars} genes")
    return adata


def run_qc(adata: sc.AnnData) -> sc.AnnData:
    """Compute QC metrics, plot them, then filter barcodes and rare genes."""
    sc.pp.filter_cells(adata, min_genes=cfg.MIN_GENES_PER_CELL)
    sc.pp.filter_genes(adata, min_cells=cfg.MIN_CELLS_PER_GENE)

    adata.var["mt"] = adata.var_names.str.startswith(cfg.MITO_PREFIX)
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
    )

    # distributions before the cut
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


def normalize_and_select(adata: sc.AnnData) -> sc.AnnData:
    """Normalise, log, pick HVGs, regress depth/mito, scale, then PCA."""
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

    # keep the full log-normalised matrix around for plotting and DE
    adata.raw = adata
    adata = adata[:, adata.var.highly_variable].copy()

    sc.pp.regress_out(adata, cfg.REGRESS_OUT)
    sc.pp.scale(adata, max_value=cfg.SCALE_MAX_VALUE)

    sc.tl.pca(adata, svd_solver="arpack", random_state=cfg.RANDOM_SEED)
    sc.pl.pca_variance_ratio(
        adata, n_pcs=50, log=True, show=False, save="_pca_variance.png"
    )
    return adata


def reduce_and_cluster(adata: sc.AnnData) -> sc.AnnData:
    """kNN graph in PCA space, UMAP for the picture, Leiden for the groups."""
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


def find_markers(adata: sc.AnnData, n_top: int = 25) -> pd.DataFrame:
    """One-vs-rest Wilcoxon test per cluster; writes results/marker_genes.csv."""
    sc.tl.rank_genes_groups(
        adata, groupby="leiden", method="wilcoxon",
        use_raw=True, n_genes=n_top,
    )
    sc.pl.rank_genes_groups(
        adata, n_genes=20, sharey=False, show=False, save="_marker_genes.png"
    )

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


def annotate_clusters(adata: sc.AnnData) -> sc.AnnData:
    """Score every cluster against the marker sets in config, take the best."""
    raw = adata.raw.to_adata()
    available = {g for g in raw.var_names}
    clusters = sorted(adata.obs["leiden"].cat.categories, key=int)

    all_markers = sorted({g for genes in cfg.MARKER_SETS.values()
                          for g in genes if g in available})
    expr = pd.DataFrame(index=clusters, columns=all_markers, dtype=float)
    for cl in clusters:
        mask = (adata.obs["leiden"] == cl).values
        sub = raw[mask, all_markers].X
        sub = np.asarray(sub.mean(axis=0)).ravel()
        expr.loc[cl] = sub

    # z-score across clusters so "high" is relative, then average per set
    zexpr = (expr - expr.mean(axis=0)) / (expr.std(axis=0) + 1e-9)

    score = pd.DataFrame(index=clusters, columns=list(cfg.MARKER_SETS), dtype=float)
    for ctype, genes in cfg.MARKER_SETS.items():
        present = [g for g in genes if g in available]
        score[ctype] = zexpr[present].mean(axis=1)

    mapping = score.idxmax(axis=1).to_dict()
    best_score = score.max(axis=1).to_dict()

    mapping_df = (
        score.assign(assigned=pd.Series(mapping), confidence=pd.Series(best_score))
        .rename_axis("cluster").reset_index()
    )
    mapping_df.to_csv(cfg.MAPPING_CSV, index=False)
    print("[annotate] cluster -> cell type:")
    for cl in clusters:
        print(f"           cluster {cl:>2} -> {mapping[cl]:<18} "
              f"(score {best_score[cl]:+.2f})")

    adata.obs["cell_type"] = (
        adata.obs["leiden"].map(mapping).astype("category")
    )

    dot_genes = {ct: [g for g in genes if g in available]
                 for ct, genes in cfg.MARKER_SETS.items()}
    sc.pl.dotplot(
        adata, dot_genes, groupby="leiden", use_raw=True,
        standard_scale="var", show=False, save="_marker_dotplot.png",
    )
    sc.pl.umap(
        adata, color="cell_type", legend_loc="right margin",
        title="PBMC 3k - annotated cell types", frameon=False,
        show=False, save="_celltypes.png",
    )
    spotlight = [g for g in cfg.SPOTLIGHT_MARKERS if g in available]
    sc.pl.umap(
        adata, color=spotlight, use_raw=True, frameon=False,
        ncols=3, show=False, save="_marker_featureplots.png",
    )

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


def tidy_figures() -> None:
    """Rename Scanpy's auto-named PNGs to the sequence in config.FIGURE_RENAMES."""
    for src_name, dst_name in cfg.FIGURE_RENAMES.items():
        src = cfg.FIGURES_DIR / src_name
        if src.exists():
            src.replace(cfg.FIGURES_DIR / dst_name)
    print(f"[tidy] renamed figures -> {len(cfg.FIGURE_RENAMES)} clean filenames")
