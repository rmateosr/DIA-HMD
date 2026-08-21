#!/bin/bash
# ABOUTME: Part of the DIA-NN Level 1 pipeline toolchain.
# ABOUTME: Post-processing: peptide FASTA conversion, canonical filtering, the detection gate,
# ABOUTME: then the TP/FP classification and R analysis jobs.
set -o errexit
set -o nounset

SCRIPT_DIR="${SGE_O_WORKDIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")" && pwd)}}"
source "$SCRIPT_DIR/config.sh"

ensure_tool python3 "python/3.12.0"

# All three come from the search, not from StrictFilter -- this stage depends on that one for
# ordering only, so the search is the stage to name if any of them is absent.
require_inputs DIANN \
  Reports/report_peptidoforms.pr_matrix.tsv \
  Reports/report_peptidoforms.parquet \
  Library/library_FROM_peptidoform.parquet

# Convert DIA-NN peptidoform matrix to FASTA for canonical filtering
# Header format: >{Protein.Group}_{Stripped.Sequence}_{Precursor.Charge}
awk -F'\t' '
  NR==1 { for(i=1;i<=NF;i++) { if($i=="Protein.Group") pg=i; if($i=="Stripped.Sequence") ss=i; if($i=="Precursor.Charge") pc=i } next }
  { print ">"$pg"_"$ss"_"$pc"\n"$ss }
' Reports/report_peptidoforms.pr_matrix.tsv > peptide.fasta

./filter_canonical_peptides.sh "$PROTEOME_FILE"

# Apply the detection gate now that the variant peptide list exists: per-run Q.Value, fragment
# geometry, and the replicate requirement if one is configured. Produces pr_matrix.strict.tsv,
# which is what the R analysis and any figure script consume -- so a plot cannot show a point the
# results table rejected.
gate_args=()
[[ -n "${SAMPLE_MAP:-}" ]] && gate_args+=(-s "$SAMPLE_MAP")
python3 "$SCRIPT_DIR/gate_variant_cells.py" \
  -q "$QVALUE_GATE" \
  --fragment-min "$FRAGMENT_MIN" \
  --min-replicates "$MIN_REPLICATES" \
  --proteome "$FASTA_FILE" \
  "${gate_args[@]+"${gate_args[@]}"}"

mkdir -p Peptidomics_Results

# TP/FP classification needs a ground truth; without one there is nothing to classify against.
if [[ -n "${TRUTH_FILE:-}" ]]; then
  submit_job Classify 4 8G "" classify_job.sh
else
  echo "TRUTH_FILE is not set in config.sh -- skipping the classification stage."
fi

submit_job RHotspot 1 8G "" submit_hotspot_analysis.sh
