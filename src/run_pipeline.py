"""Runs the full PBMC 3k analysis end to end: python -m src.run_pipeline"""

import time

from . import config as cfg
from . import pipeline as pl


def main() -> None:
    t0 = time.time()

    pl.configure_scanpy()

    adata = pl.load_data()
    adata = pl.run_qc(adata)
    adata = pl.normalize_and_select(adata)
    adata = pl.reduce_and_cluster(adata)
    pl.find_markers(adata)
    adata = pl.annotate_clusters(adata)
    pl.tidy_figures()

    adata.write_h5ad(cfg.PROCESSED_H5AD)
    print(f"\n[save] wrote {cfg.PROCESSED_H5AD}")

    dt = time.time() - t0
    print(f"Done in {dt:0.1f}s. "
          f"{adata.n_obs} cells, {adata.obs['cell_type'].nunique()} cell types.")
    print(f"Figures -> {cfg.FIGURES_DIR}")
    print(f"Tables  -> {cfg.RESULTS_DIR}")


if __name__ == "__main__":
    main()
