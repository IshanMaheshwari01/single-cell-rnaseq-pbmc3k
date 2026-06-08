"""
config.py
=========
Central configuration for the PBMC 3k single-cell RNA-seq pipeline.

Keeping every "magic number" in one place is a deliberate design choice:
in real single-cell work you will re-run the pipeline many times while you
tune the quality-control thresholds, the number of principal components, and
the clustering resolution. Having them here (instead of scattered through the
code) makes the analysis reproducible and easy to audit.
"""

from pathlib import Path

# --------------------------------------------------------------------------- #
# Paths
# --------------------------------------------------------------------------- #
# Resolve paths relative to the project root so the pipeline works no matter
# what directory you launch it from.
PROJECT_ROOT = Path(__file__).resolve().parents[1]

DATA_DIR = PROJECT_ROOT / "data" / "pbmc3k" / "filtered_gene_bc_matrices" / "hg19"
FIGURES_DIR = PROJECT_ROOT / "figures"
RESULTS_DIR = PROJECT_ROOT / "results"

PROCESSED_H5AD = RESULTS_DIR / "pbmc3k_processed.h5ad"
MARKER_CSV = RESULTS_DIR / "marker_genes.csv"
MAPPING_CSV = RESULTS_DIR / "cluster_celltype_mapping.csv"
COUNTS_CSV = RESULTS_DIR / "celltype_counts.csv"

# Make sure output folders exist.
FIGURES_DIR.mkdir(parents=True, exist_ok=True)
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# --------------------------------------------------------------------------- #
# Reproducibility
# --------------------------------------------------------------------------- #
# A single seed used everywhere (PCA, neighbours, UMAP, Leiden). Stochastic
# steps such as UMAP and community detection will give slightly different
# pictures on every run unless the seed is fixed.
RANDOM_SEED = 0

# --------------------------------------------------------------------------- #
# Quality-control thresholds
# --------------------------------------------------------------------------- #
# A "cell" in droplet-based scRNA-seq is really a barcode. Some barcodes are
# empty droplets (almost no genes), some are doublets (two cells in one
# droplet -> abnormally many genes), and some are dying cells (high fraction of
# mitochondrial transcripts leaking out). These thresholds remove those.
MIN_GENES_PER_CELL = 200   # drop barcodes expressing < 200 genes (empty / debris)
MIN_CELLS_PER_GENE = 3     # drop genes seen in < 3 cells (noise, not informative)
MAX_GENES_PER_CELL = 2500  # drop barcodes expressing > 2500 genes (likely doublets)
MAX_PCT_MITO = 5.0         # drop cells with > 5% mitochondrial counts (stressed/dying)
MITO_PREFIX = "MT-"        # human mitochondrial gene symbols start with "MT-"

# --------------------------------------------------------------------------- #
# Normalisation
# --------------------------------------------------------------------------- #
# Library-size normalisation: scale every cell to the same total counts so that
# a cell sequenced more deeply does not look "more expressed" than its neighbour.
TARGET_SUM = 1e4  # counts-per-10k (a.k.a. CP10K)

# --------------------------------------------------------------------------- #
# Highly variable genes (feature selection)
# --------------------------------------------------------------------------- #
# Most of the ~14k detected genes are housekeeping / uninformative. We keep the
# genes whose variability across cells exceeds what we would expect from
# technical noise alone. Seurat's "dispersion" recipe is used here.
HVG_MIN_MEAN = 0.0125
HVG_MAX_MEAN = 3.0
HVG_MIN_DISP = 0.5

# Technical covariates to regress out before scaling (removes the effect of
# sequencing depth and mitochondrial fraction on downstream structure).
REGRESS_OUT = ["total_counts", "pct_counts_mt"]
SCALE_MAX_VALUE = 10.0  # clip scaled expression to +/- 10 SD to limit outliers

# --------------------------------------------------------------------------- #
# Dimensionality reduction & clustering
# --------------------------------------------------------------------------- #
N_PCS = 40          # number of principal components carried into the graph
N_NEIGHBORS = 10    # k for the k-nearest-neighbour graph
LEIDEN_RESOLUTION = 1.0  # higher -> more, smaller clusters

# --------------------------------------------------------------------------- #
# Cell-type annotation
# --------------------------------------------------------------------------- #
# Canonical PBMC marker genes. These are textbook lineage markers for the major
# immune populations you expect to find in peripheral blood. The annotation step
# scores each cluster against every set and assigns the best match.
MARKER_SETS = {
    "CD4 T cells":       ["IL7R", "CD3D", "CD3E", "CCR7", "LDHB"],
    "CD8 T cells":       ["CD8A", "CD8B", "CD3D", "CD3E", "GZMK"],
    "NK cells":          ["GNLY", "NKG7", "KLRD1", "NCAM1", "GZMB"],
    "B cells":           ["MS4A1", "CD79A", "CD79B", "CD19"],
    "CD14+ Monocytes":   ["CD14", "LYZ", "S100A8", "S100A9"],
    "FCGR3A+ Monocytes": ["FCGR3A", "MS4A7", "CDKN1C"],
    "Dendritic cells":   ["FCER1A", "CST3", "CLEC10A"],
    "Megakaryocytes":    ["PPBP", "PF4", "ITGA2B"],
}

# A compact set of one-marker-per-lineage genes used for the UMAP feature plots.
SPOTLIGHT_MARKERS = [
    "CD3D",    # pan T cell
    "IL7R",    # CD4 T
    "CD8A",    # CD8 T
    "NKG7",    # NK
    "MS4A1",   # B
    "CD14",    # CD14 monocyte
    "FCGR3A",  # FCGR3A monocyte
    "FCER1A",  # dendritic
    "PPBP",    # megakaryocyte / platelet
]

# --------------------------------------------------------------------------- #
# Figure file names
# --------------------------------------------------------------------------- #
# Scanpy prepends the plot type to every file it saves (e.g. "violin_...").
# We rename them to a clean, numbered sequence after the run so the figures/
# folder reads top-to-bottom like the analysis itself. Map: scanpy name -> clean.
FIGURE_RENAMES = {
    "violin_qc_before.png":                    "01_qc_violin_before.png",
    "scatter_qc_counts_vs_mito.png":           "02_qc_scatter_counts_vs_mito.png",
    "scatter_qc_counts_vs_genes.png":          "03_qc_scatter_counts_vs_genes.png",
    "violin_qc_after.png":                     "04_qc_violin_after.png",
    "filter_genes_dispersion_hvg.png":         "05_highly_variable_genes.png",
    "pca_variance_ratio_pca_variance.png":     "06_pca_variance_ratio.png",
    "umap_clusters.png":                       "07_umap_leiden_clusters.png",
    "umap_qc_on_umap.png":                     "08_umap_qc_metrics.png",
    "rank_genes_groups_leiden_marker_genes.png": "09_rank_genes_groups.png",
    "dotplot__marker_dotplot.png":             "10_marker_dotplot.png",
    "umap_celltypes.png":                      "11_umap_celltypes.png",
    "umap_marker_featureplots.png":            "12_umap_marker_featureplots.png",
    "celltype_proportions.png":                "13_celltype_proportions.png",
}
