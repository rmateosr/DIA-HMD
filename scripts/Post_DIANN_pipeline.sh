#!/bin/bash
# ABOUTME: Part of the DIA-NN Level 1 pipeline toolchain.
# ABOUTME: Post-processing: peptide FASTA conversion, canonical filtering, then R analysis jobs.
set -o errexit
set -o nounset

SCRIPT_DIR="${SGE_O_WORKDIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")" && pwd)}}"
source "$SCRIPT_DIR/config.sh"

# Convert DIA-NN peptidoform matrix to FASTA for canonical filtering
# Header format: >{Protein.Group}_{Stripped.Sequence}_{Precursor.Charge}
awk -F'\t' '
  NR==1 { for(i=1;i<=NF;i++) { if($i=="Protein.Group") pg=i; if($i=="Stripped.Sequence") ss=i; if($i=="Precursor.Charge") pc=i } next }
  { print ">"$pg"_"$ss"_"$pc"\n"$ss }
' Reports/report_peptidoforms.pr_matrix.strict.tsv > peptide.fasta

./filter_canonical_peptides.sh "$PROTEOME_FILE"

mkdir -p Peptidomics_Results
submit_job RHotspot 1 8G "" submit_hotspot_analysis.sh
