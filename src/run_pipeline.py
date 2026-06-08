"""
run_pipeline.py
===============
End-to-end driver for the PBMC 3k single-cell RNA-seq analysis.

Run it from the project root with:

    python -m src.run_pipeline

It executes every stage in order, writes all figures to ./figures, all tables
to ./results, and saves the fully processed, annotated AnnData object to
./results/pbmc3k_processed.h5ad so you can re-load it instantly later:

    import scanpy as sc
    adata = sc.read_h5ad("results/pbmc3k_processed.h5ad")
"""

import time

from . import config as cfg
from . import pipeline as pl


def main() -> None:
    t0 = time.time()
    print("=" * 64)
    print("PBMC 3k single-cell RNA-seq pipeline")
    print("=" * 64)

    pl.configure_scanpy()

    adata = pl.load_data()
    adata = pl.run_qc(adata)
    adata = pl.normalize_and_select(adata)
    adata = pl.reduce_and_cluster(adata)
    pl.find_markers(adata)
    adata = pl.annotate_clusters(adata)
    pl.tidy_figures()

    # Persist the final object for instant re-loading / sharing.
    adata.write_h5ad(cfg.PROCESSED_H5AD)
    print(f"\n[save] wrote {cfg.PROCESSED_H5AD}")

    dt = time.time() - t0
    print("=" * 64)
    print(f"Done in {dt:0.1f}s. "
          f"{adata.n_obs} cells, {adata.obs['cell_type'].nunique()} cell types.")
    print(f"Figures -> {cfg.FIGURES_DIR}")
    print(f"Tables  -> {cfg.RESULTS_DIR}")
    print("=" * 64)


if __name__ == "__main__":
    main()
