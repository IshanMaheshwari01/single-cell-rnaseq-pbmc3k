"""All thresholds and paths for the PBMC 3k pipeline in one place."""

from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]

DATA_DIR = PROJECT_ROOT / "data" / "pbmc3k" / "filtered_gene_bc_matrices" / "hg19"
FIGURES_DIR = PROJECT_ROOT / "figures"
RESULTS_DIR = PROJECT_ROOT / "results"

PROCESSED_H5AD = RESULTS_DIR / "pbmc3k_processed.h5ad"
MARKER_CSV = RESULTS_DIR / "marker_genes.csv"
MAPPING_CSV = RESULTS_DIR / "cluster_celltype_mapping.csv"
COUNTS_CSV = RESULTS_DIR / "celltype_counts.csv"

FIGURES_DIR.mkdir(parents=True, exist_ok=True)
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# one seed everywhere (PCA, neighbours, UMAP, Leiden)
RANDOM_SEED = 0

# QC thresholds
MIN_GENES_PER_CELL = 200   # empty droplets / debris
MIN_CELLS_PER_GENE = 3     # genes seen in too few cells to be informative
MAX_GENES_PER_CELL = 2500  # likely doublets
MAX_PCT_MITO = 5.0         # stressed or dying cells
MITO_PREFIX = "MT-"

# normalisation
TARGET_SUM = 1e4  # counts per 10k

# highly variable genes (Seurat dispersion recipe)
HVG_MIN_MEAN = 0.0125
HVG_MAX_MEAN = 3.0
HVG_MIN_DISP = 0.5

# technical covariates removed before scaling
REGRESS_OUT = ["total_counts", "pct_counts_mt"]
SCALE_MAX_VALUE = 10.0

# dimensionality reduction & clustering
N_PCS = 40
N_NEIGHBORS = 10
LEIDEN_RESOLUTION = 1.0

# canonical PBMC markers, used to name the clusters
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

# one marker per lineage, for the UMAP feature plots
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

# Scanpy names saved plots after the plot type; rename them afterwards so
# figures/ reads in analysis order. Map: scanpy name -> clean name.
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
