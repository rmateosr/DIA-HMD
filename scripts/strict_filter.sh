#!/bin/bash
# ABOUTME: Summarises the DIA-NN report as a per-peptide table of best q-values and runs detected.
# ABOUTME: Reporting only -- the per-sample gate is applied later, by gate_variant_cells.py.
set -euo pipefail

SCRIPT_DIR="${SGE_O_WORKDIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")" && pwd)}}"
source "$SCRIPT_DIR/config.sh"

ensure_tool python3 "python/3.12.0"

echo "=== Precursor summary ==="

python3 "$SCRIPT_DIR/extract_strict_precursors.py" \
  Reports/report_peptidoforms.parquet \
  -o Reports/strict_precursors_peptide_list.tsv

echo "=== Precursor summary complete ==="
