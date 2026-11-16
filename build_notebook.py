"""
Rebuilds notebooks/pbmc3k_analysis.ipynb from the markdown + code cells below,
so the notebook stays reviewable as plain Python.

    python build_notebook.py
    jupyter nbconvert --to notebook --execute --inplace \
        --ExecutePreprocessor.timeout=900 notebooks/pbmc3k_analysis.ipynb
"""

import nbformat as nbf
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell

nb = new_notebook()
cells = []


def md(text):
    cells.append(new_markdown_cell(text.strip("\n")))


def code(text):
    cells.append(new_code_cell(text.strip("\n")))


md(r"""
# Single-cell RNA-seq of 2,700 PBMCs

**Author:** [IshanMaheshwari01](https://github.com/IshanMaheshwari01)
**Tools:** Python · [Scanpy](https://scanpy.readthedocs.io) · AnnData · Leiden · UMAP

Bulk RNA-seq averages expression over a whole tissue. Single-cell RNA-seq measures
each cell separately, so we can work out which cell types are in a sample and
study them on their own. The catch is the data is sparse and noisy.

The dataset is the classic **PBMC 3k** set from 10x Genomics: 2,700 peripheral
blood mononuclear cells from one healthy donor. PBMCs are the immune cells in
blood, so we know roughly which cell types to expect.

The goal: take the raw gene-by-cell count matrix and, with no prior labels,
recover the known immune cell types and the marker genes that define them.

1. QC: remove empty droplets, doublets and dying cells
2. Normalisation
3. Highly variable genes
4. PCA, then UMAP
5. Leiden clustering
6. Marker genes
7. Annotation
""")

md(r"""
## 0. Setup

Fixed `random_state` everywhere, otherwise UMAP and clustering come out
slightly different on every run.
""")

code(r"""
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from pathlib import Path

SEED = 0
np.random.seed(SEED)

sc.settings.verbosity = 1
sc.settings.set_figure_params(dpi=90, facecolor="white", frameon=False)

# data lives at <project_root>/data/... ; this notebook lives in notebooks/
DATA_DIR = Path("..") / "data" / "pbmc3k" / "filtered_gene_bc_matrices" / "hg19"
print("scanpy", sc.__version__)
""")

md(r"""
## 1. Load the data

10x Cell Ranger outputs three files: `matrix.mtx` (the counts), `genes.tsv`
(gene names, rows) and `barcodes.tsv` (cell barcodes, columns). Scanpy reads
them into an AnnData object, which keeps the expression matrix (`.X`), cell
metadata (`.obs`) and gene metadata (`.var`) together.
""")

code(r"""
adata = sc.read_10x_mtx(DATA_DIR, var_names="gene_symbols", cache=True)
adata.var_names_make_unique()
adata
""")

md(r"""
2,700 cells and 32,738 genes. The matrix is ~97% zeros, which is normal and
is why the steps below are needed.
""")

md(r"""
## 2. Quality control (QC)

A "cell" here is a droplet barcode, and not every barcode is a healthy single
cell. Three diagnostics per cell:

- `n_genes_by_counts`: too few means an empty droplet, too many a doublet
- `total_counts`: sequencing depth
- `pct_counts_mt`: mitochondrial fraction. High values mean a stressed or
  dying cell, since cytoplasmic mRNA leaks out and the mitochondrial fraction
  goes up.
""")

code(r"""
# light pre-filter: drop near-empty barcodes and genes seen in <3 cells
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

# flag mitochondrial genes and compute standard QC metrics
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None,
                           log1p=False, inplace=True)

sc.pl.violin(adata, ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
             jitter=0.4, multi_panel=True)
""")

md(r"""
Cells with very high mitochondrial percentage or an unusually large number of
genes are the ones to exclude.
""")

code(r"""
sc.pl.scatter(adata, x="total_counts", y="pct_counts_mt")
sc.pl.scatter(adata, x="total_counts", y="n_genes_by_counts")
""")

code(r"""
n_before = adata.n_obs
adata = adata[adata.obs.n_genes_by_counts < 2500, :].copy()  # remove doublets
adata = adata[adata.obs.pct_counts_mt < 5, :].copy()         # remove dying cells
print(f"kept {adata.n_obs} / {n_before} cells after QC")
""")

md(r"""
## 3. Normalisation & log-transform

A cell sequenced twice as deeply shows roughly twice the counts everywhere.
That's technical, not biological. Two steps: scale every cell to the same
total (10k), then `log(1+x)` so a handful of highly expressed genes don't
dominate distances between cells.
""")

code(r"""
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
""")

md(r"""
## 4. Highly variable genes

Most of the remaining ~13,700 genes are housekeeping genes that look the same
in every cell. Keep only the highly variable ones, freeze the full
log-normalised matrix in `adata.raw` for plotting later, and work with the
HVGs from here on.
""")

code(r"""
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
print(f"{int(adata.var.highly_variable.sum())} highly variable genes selected")
sc.pl.highly_variable_genes(adata)
""")

code(r"""
adata.raw = adata                                  # keep all genes for plotting
adata = adata[:, adata.var.highly_variable].copy() # work with HVGs only
""")

md(r"""
## 5. Scaling & PCA

Regress out total counts and mitochondrial fraction, scale each gene to mean 0
/ unit variance (clipped at ±10), then PCA. The elbow in the variance plot
says how many PCs carry real signal.
""")

code(r"""
sc.pp.regress_out(adata, ["total_counts", "pct_counts_mt"])
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver="arpack", random_state=SEED)
sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True)
""")

md(r"""
## 6. Neighbourhood graph, UMAP & clustering

kNN graph on the top 40 PCs, a UMAP embedding of that graph for looking at,
and Leiden clustering on the graph itself (not the 2-D coordinates).
""")

code(r"""
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40, random_state=SEED)
sc.tl.umap(adata, random_state=SEED)
sc.tl.leiden(adata, resolution=1.0, random_state=SEED,
             flavor="igraph", n_iterations=2, directed=False)
print(f"Leiden found {adata.obs['leiden'].nunique()} clusters")
sc.pl.umap(adata, color="leiden", legend_loc="on data", title="Leiden clusters")
""")

md(r"""
Colouring the same UMAP by the QC metrics is a sanity check: if a cluster were
defined by sequencing depth or mitochondrial content rather than biology, it
would light up here.
""")

code(r"""
sc.pl.umap(adata, color=["n_genes_by_counts", "total_counts", "pct_counts_mt"])
""")

md(r"""
## 7. Marker genes

One-vs-rest differential expression (Wilcoxon rank-sum) per cluster, to find
the genes most specific to each one.
""")

code(r"""
sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon", use_raw=True)
sc.pl.rank_genes_groups(adata, n_genes=20, sharey=False)
""")

code(r"""
# top 5 markers per cluster, as a tidy table
top = sc.get.rank_genes_groups_df(adata, group=None)
top.groupby("group").head(5).reset_index(drop=True)
""")

md(r"""
## 8. Do the markers match known immune cell types?

Dotplot of canonical PBMC markers: colour is mean expression, dot size is the
fraction of cells expressing the gene. A clean diagonal means each cluster has
its own marker signature.

| Cell type | Canonical markers |
|-----------|-------------------|
| CD4 T cells | IL7R, CD3D, CCR7 |
| CD8 T cells | CD8A, CD8B, GZMK |
| NK cells | GNLY, NKG7, KLRD1 |
| B cells | MS4A1, CD79A, CD79B |
| CD14+ Monocytes | CD14, LYZ, S100A8/9 |
| FCGR3A+ Monocytes | FCGR3A, MS4A7 |
| Dendritic cells | FCER1A, CST3 |
| Megakaryocytes | PPBP, PF4 |
""")

code(r"""
marker_sets = {
    "CD4 T":        ["IL7R", "CD3D", "CD3E", "CCR7", "LDHB"],
    "CD8 T":        ["CD8A", "CD8B", "CD3D", "CD3E", "GZMK"],
    "NK":           ["GNLY", "NKG7", "KLRD1", "NCAM1", "GZMB"],
    "B":            ["MS4A1", "CD79A", "CD79B", "CD19"],
    "CD14+ Mono":   ["CD14", "LYZ", "S100A8", "S100A9"],
    "FCGR3A+ Mono": ["FCGR3A", "MS4A7", "CDKN1C"],
    "DC":           ["FCER1A", "CST3", "CLEC10A"],
    "Megakaryo":    ["PPBP", "PF4", "ITGA2B"],
}
sc.pl.dotplot(adata, marker_sets, groupby="leiden",
              use_raw=True, standard_scale="var")
""")

md(r"""
## 9. Annotate the clusters

Score each cluster against the marker sets (mean z-scored expression) and take
the best match, rather than hard-coding "cluster 3 = CD8 T", so the labels
survive cluster renumbering between software versions.
""")

code(r"""
raw = adata.raw.to_adata()
clusters = sorted(adata.obs["leiden"].cat.categories, key=int)
genes = sorted({g for v in marker_sets.values() for g in v if g in set(raw.var_names)})

# mean expression of each marker gene per cluster, z-scored across clusters
expr = pd.DataFrame(
    {cl: np.asarray(raw[(adata.obs.leiden == cl).values, genes].X.mean(0)).ravel()
     for cl in clusters},
    index=genes,
).T
z = (expr - expr.mean(0)) / (expr.std(0) + 1e-9)

score = pd.DataFrame(
    {ct: z[[g for g in gs if g in genes]].mean(1) for ct, gs in marker_sets.items()}
)
mapping = score.idxmax(1).to_dict()
mapping
""")

code(r"""
full_names = {
    "CD4 T": "CD4 T cells", "CD8 T": "CD8 T cells", "NK": "NK cells",
    "B": "B cells", "CD14+ Mono": "CD14+ Monocytes",
    "FCGR3A+ Mono": "FCGR3A+ Monocytes", "DC": "Dendritic cells",
    "Megakaryo": "Megakaryocytes",
}
adata.obs["cell_type"] = (
    adata.obs["leiden"].map(lambda c: full_names[mapping[c]]).astype("category")
)
sc.pl.umap(adata, color="cell_type", legend_loc="right margin",
           title="PBMC 3k: annotated cell types")
""")

md(r"""
## 10. Confirm with single-marker feature plots

One canonical marker per lineage on the UMAP; each should light up exactly
one region.
""")

code(r"""
sc.pl.umap(adata, color=["CD3D", "IL7R", "CD8A", "NKG7", "MS4A1",
                         "CD14", "FCGR3A", "FCER1A", "PPBP"],
           use_raw=True, ncols=3)
""")

md(r"""
## 11. Cell-type composition
""")

code(r"""
counts = adata.obs["cell_type"].value_counts()
ax = counts.sort_values().plot(kind="barh", figsize=(7, 4), color="#4C72B0")
ax.set_xlabel("number of cells"); ax.set_title("Cell-type composition of PBMC 3k")
plt.tight_layout(); plt.show()
counts.to_frame("n_cells")
""")

md(r"""
## Conclusion

From a raw, unlabelled 2,700 x 32,738 count matrix the pipeline recovered
eight immune populations: CD4 and CD8 T cells, NK cells, B cells, two monocyte
subsets, dendritic cells and megakaryocytes. T cells dominate, as in real
blood, and even the 13-cell megakaryocyte population came out cleanly.

Obvious next steps: cell-cycle regression, doublet detection (Scrublet),
pseudotime, batch integration (Harmony / scVI) and reference-based annotation
(CellTypist).
""")

nb["cells"] = cells
nb["metadata"] = {
    "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
    "language_info": {"name": "python"},
}

out = "notebooks/pbmc3k_analysis.ipynb"
with open(out, "w") as f:
    nbf.write(nb, f)
print(f"wrote {out} with {len(cells)} cells")
