#!/usr/bin/env bash
# ===========================================================================
# download_data.sh  —  fetch the 10x Genomics PBMC 3k dataset.
#
# This script reconstructs the standard Cell Ranger v2 layout that the pipeline
# expects:
#
#     data/pbmc3k/filtered_gene_bc_matrices/hg19/
#         ├── matrix.mtx
#         ├── genes.tsv
#         └── barcodes.tsv
#
# It tries the official 10x Genomics download first. If that host is
# unreachable (firewalls, mirrors going away, etc.) it falls back to a GitHub
# mirror that stores the same unmodified triple, and rebuilds the layout.
#
# Usage:
#     bash data/download_data.sh
# ===========================================================================
set -euo pipefail

# Resolve paths relative to this script so it works from any directory.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
DATA_DIR="$PROJECT_ROOT/data/pbmc3k"
HG19_DIR="$DATA_DIR/filtered_gene_bc_matrices/hg19"

OFFICIAL_URL="https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz"
MIRROR_URL="https://codeload.github.com/Amartya101/PBMC3k_Data/tar.gz/refs/heads/main"

# --------------------------------------------------------------------------- #
# 0. Skip if the data is already in place.
# --------------------------------------------------------------------------- #
if [[ -f "$HG19_DIR/matrix.mtx" && -f "$HG19_DIR/genes.tsv" && -f "$HG19_DIR/barcodes.tsv" ]]; then
    echo "[download_data] Data already present at:"
    echo "                $HG19_DIR"
    echo "[download_data] Nothing to do. (Delete data/pbmc3k/ to force a re-download.)"
    exit 0
fi

mkdir -p "$DATA_DIR"
TMP_DIR="$(mktemp -d)"
trap 'rm -rf "$TMP_DIR"' EXIT

verify() {
    # Sanity-check the reconstructed layout: ~2,700 barcodes, ~32,738 genes.
    local n_bc n_genes
    n_bc=$(wc -l < "$HG19_DIR/barcodes.tsv" | tr -d ' ')
    n_genes=$(wc -l < "$HG19_DIR/genes.tsv" | tr -d ' ')
    echo "[download_data] Verifying: $n_bc barcodes, $n_genes genes."
    if [[ "$n_bc" -lt 2000 || "$n_genes" -lt 20000 ]]; then
        echo "[download_data] WARNING: counts look smaller than expected." >&2
    fi
}

# --------------------------------------------------------------------------- #
# 1. Try the official 10x Genomics download.
#    The tarball extracts to filtered_gene_bc_matrices/hg19/{matrix.mtx,
#    genes.tsv, barcodes.tsv}, i.e. exactly our target layout.
# --------------------------------------------------------------------------- #
echo "[download_data] Attempting official 10x Genomics download..."
if curl -fL --connect-timeout 20 -o "$TMP_DIR/pbmc3k.tar.gz" "$OFFICIAL_URL"; then
    echo "[download_data] Downloaded from 10x. Extracting..."
    tar -xzf "$TMP_DIR/pbmc3k.tar.gz" -C "$DATA_DIR"
    if [[ -f "$HG19_DIR/matrix.mtx" && -f "$HG19_DIR/genes.tsv" && -f "$HG19_DIR/barcodes.tsv" ]]; then
        echo "[download_data] Success (official 10x source)."
        verify
        exit 0
    fi
    echo "[download_data] Official archive had an unexpected layout; falling back." >&2
else
    echo "[download_data] Official host unreachable; falling back to GitHub mirror." >&2
fi

# --------------------------------------------------------------------------- #
# 2. Fallback: GitHub mirror (Amartya101/PBMC3k_Data).
#    Stores the same unmodified triple as:
#        PBMC3k_Data-main/PBMC3k/10X_PBMC3k_barcodes.tsv
#        PBMC3k_Data-main/PBMC3k/10X_PBMC3k_features.tsv
#        PBMC3k_Data-main/PBMC3k/10X_PBMC3k_matrix.mtx.gz
#    We rebuild the standard hg19 layout from it.
# --------------------------------------------------------------------------- #
echo "[download_data] Downloading from GitHub mirror..."
curl -fL --connect-timeout 30 -o "$TMP_DIR/mirror.tar.gz" "$MIRROR_URL"
echo "[download_data] Extracting mirror..."
tar -xzf "$TMP_DIR/mirror.tar.gz" -C "$TMP_DIR"

SRC="$TMP_DIR/PBMC3k_Data-main/PBMC3k"
if [[ ! -d "$SRC" ]]; then
    echo "[download_data] ERROR: unexpected mirror layout; aborting." >&2
    exit 1
fi

mkdir -p "$HG19_DIR"
gunzip -c "$SRC/10X_PBMC3k_matrix.mtx.gz" > "$HG19_DIR/matrix.mtx"
cp "$SRC/10X_PBMC3k_features.tsv" "$HG19_DIR/genes.tsv"
cp "$SRC/10X_PBMC3k_barcodes.tsv" "$HG19_DIR/barcodes.tsv"

echo "[download_data] Success (GitHub mirror)."
verify
echo "[download_data] Data ready at: $HG19_DIR"
