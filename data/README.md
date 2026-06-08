# Data — PBMC 3k (10x Genomics)

This folder holds the input dataset for the analysis. The actual data files are
**not committed to git** (they are re-downloadable and a bit large); run the
download script to fetch them:

```bash
bash data/download_data.sh
```

After running it you will have the standard Cell Ranger layout:

```
data/pbmc3k/filtered_gene_bc_matrices/hg19/
├── matrix.mtx      # sparse gene × cell count matrix (Matrix Market format)
├── genes.tsv       # gene IDs + symbols (one row per gene)
└── barcodes.tsv    # cell barcodes (one row per cell)
```

## What is this dataset?

**PBMC 3k** is one of the most widely used reference datasets in single-cell
genomics. It contains **2,700 peripheral blood mononuclear cells (PBMCs)** from
a single healthy human donor, sequenced on the 10x Genomics Chromium platform
and processed with Cell Ranger against the **hg19** reference. PBMCs are the
immune cells found in blood — T cells, B cells, NK cells, monocytes, dendritic
cells — which makes the dataset ideal for learning and benchmarking: the cell
types that *should* appear are well known, so you can check whether an analysis
recovers the right biology.

Raw matrix dimensions (before any filtering): **32,738 genes × 2,700 cells**.

## Provenance

The pipeline is written against the canonical 10x Genomics release. The
download script tries two sources, in order:

1. **Official 10x Genomics** —
   `https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz`
   This is the authoritative source and the preferred path. Its tarball expands
   directly into the `filtered_gene_bc_matrices/hg19/` layout above.

2. **GitHub mirror (fallback)** —
   [`Amartya101/PBMC3k_Data`](https://github.com/Amartya101/PBMC3k_Data),
   used only if the official host is unreachable from your network. It stores
   the same unmodified barcodes / features / matrix triple; the script rebuilds
   the standard `hg19/` layout from it. The triple was verified to match the
   canonical reference exactly (same 2,700 barcodes beginning
   `AAACATACAACCAC-1`, same 32,738 genes, identical nonzero count).

The same dataset also ships inside Scanpy as `scanpy.datasets.pbmc3k()` and in
the Seurat "Guided Clustering" tutorial; all three are the same underlying 10x
release.

## Licence / terms

The PBMC 3k dataset is produced and distributed by **10x Genomics** and is
subject to **10x Genomics' terms of use**. It is provided by 10x for research
and educational purposes. This repository redistributes none of the data
itself — it only downloads it at build time from the sources above. The MIT
licence in the repository root covers the *code*, not the data.
