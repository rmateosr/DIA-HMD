#!/bin/bash
# ABOUTME: Applies strict q-value filtering to DIA-NN output before canonical peptide filtering.
# ABOUTME: Runs extract_strict_precursors.py then filter_pr_matrix.py to produce pr_matrix.strict.tsv.
set -euo pipefail

SCRIPT_DIR="${SGE_O_WORKDIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")" && pwd)}}"
source "$SCRIPT_DIR/config.sh"

ensure_tool python3 "python/3.12.0"

echo "=== Strict precursor filtering ==="

python3 "$SCRIPT_DIR/extract_strict_precursors.py" \
  Reports/report_peptidoforms.parquet \
  -o Reports/strict_precursors_peptide_list.tsv

python3 "$SCRIPT_DIR/filter_pr_matrix.py" \
  -m Reports/report_peptidoforms.pr_matrix.tsv \
  -p Reports/strict_precursors_peptide_list.tsv \
  -o Reports/report_peptidoforms.pr_matrix.strict.tsv

echo "=== Strict filtering complete ==="
