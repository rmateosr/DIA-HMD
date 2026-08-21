#!/bin/bash
# ABOUTME: Classifies variant-peptide detections as TP/FP against the cohort ground truth.
# ABOUTME: Applies the per-run Q.Value gate, the fragment-geometry filter and the replicate rule.
# Resources: 4 slots x 8G, ~1 min (measured 3.3G peak RSS / 11 s on a 1.5G report; the slots buy
# the vmem ceiling, not parallelism -- pyarrow reserves far more than it resides).
[[ -n "${DEBUG:-}" ]] && set -xv
set -euo pipefail

SCRIPT_DIR="${SGE_O_WORKDIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")" && pwd)}}"
source "$SCRIPT_DIR/config.sh"

ensure_tool python3 "python/3.12.0"

require_inputs PostDIANN non_canonical_peptide_headers.txt
require_inputs DIANN Reports/report_peptidoforms.parquet Library/library_FROM_peptidoform.parquet

mkdir -p Peptidomics_Results

echo "=== Hotspot detection classification ==="

# Optional cohort description; each flag is only passed when the file is configured, because the
# script's own default for all three is "absent", not "missing file".
extra=()
[[ -n "${SAMPLE_MAP:-}"   ]] && extra+=(-s "$SAMPLE_MAP")
[[ -n "${ALIASES_FILE:-}" ]] && extra+=(--aliases "$ALIASES_FILE")
[[ -n "${POOLS_FILE:-}"   ]] && extra+=(--pools "$POOLS_FILE")

# -p is the FASTA DIA-NN actually searched, so the accession->gene map cannot drift from it.
# -t comes from config.sh: the recall denominator has to be the ground truth that matches the
# library that was searched.
python3 "$SCRIPT_DIR/classify_hotspot_detections.py" \
  -p "$FASTA_FILE" \
  -t "$TRUTH_FILE" \
  -q "$QVALUE_GATE" \
  --fragment-min "$FRAGMENT_MIN" \
  --min-replicates "$MIN_REPLICATES" \
  "${extra[@]+"${extra[@]}"}" \
  > Peptidomics_Results/hotspot_detection_classification_summary.txt

echo "=== Classification complete ==="
