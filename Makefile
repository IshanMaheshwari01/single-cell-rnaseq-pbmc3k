# ===========================================================================
# Makefile  —  convenience wrapper around the project's common commands.
#
# Usage:
#     make setup      install Python dependencies (pip)
#     make data       download the 10x PBMC 3k dataset
#     make run        run the end-to-end pipeline (writes figures/ + results/)
#     make notebook   execute the analysis notebook in place
#     make all        data + run
#     make clean      remove generated outputs and caches
#     make help       list available targets
#
# Tip: run `make setup` inside a fresh virtual environment, e.g.
#     python3 -m venv .venv && source .venv/bin/activate && make setup
# ===========================================================================

PYTHON ?= python3
NOTEBOOK = notebooks/pbmc3k_analysis.ipynb

.PHONY: help setup data run notebook all clean

help:
	@echo "Available targets:"
	@echo "  setup     - pip install -r requirements.txt"
	@echo "  data      - download the 10x PBMC 3k dataset"
	@echo "  run       - run the full pipeline (figures/ + results/)"
	@echo "  notebook  - execute the analysis notebook in place"
	@echo "  all       - data + run"
	@echo "  clean     - remove generated outputs and caches"

setup:
	$(PYTHON) -m pip install -r requirements.txt

data:
	bash data/download_data.sh

run:
	$(PYTHON) -m src.run_pipeline

notebook:
	jupyter nbconvert --to notebook --execute --inplace $(NOTEBOOK)

all: data run

clean:
	rm -rf cache notebooks/cache **/__pycache__ src/__pycache__ .ipynb_checkpoints
	rm -f results/*.h5ad
	@echo "Cleaned caches and large processed objects."
	@echo "Note: figures/ and results/*.csv are kept (they are committed)."
