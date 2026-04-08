#!/bin/bash
# ABOUTME: Part of the DIA-NN Level 1 pipeline toolchain.
# ABOUTME: Filters peptide.fasta against the canonical proteome; writes non-canonical headers.
start=$(date +%s)

[[ -n "${DEBUG:-}" ]] && set -xv
set -o errexit
set -o nounset

QUERY="peptide.fasta"
DB="$1"
NOT_PRESENT="non_canonical_peptide_headers.txt"

# Flatten the proteome to a single string so peptide sequences can be found with grep -F
awk '/^>/{if(NR>1) printf "X"} !/^>/{printf "%s",$0} END{printf "\n"}' "$DB" > db_seq.txt

# Convert FASTA to header<TAB>sequence TSV
awk 'BEGIN{RS=">"; ORS=""} NR>1 {n=split($0, lines, "\n"); header=lines[1]; seq=""; for (i=2; i<=n; i++) seq=seq lines[i]; print header "\t" seq "\n"}' "$QUERY" > tmp_query.tsv

> "$NOT_PRESENT"

# Split input into chunks and search in parallel using background jobs
NJOBS=8
split -n l/$NJOBS tmp_query.tsv _chunk_
pids=()
for chunk in _chunk_*; do
  (while IFS=$'\t' read -r header seq; do
    if ! grep -m 1 -qF "$seq" db_seq.txt; then
      echo "$header"
    fi
  done < "$chunk" > "${chunk}.out") &
  pids+=($!)
done
for pid in "${pids[@]}"; do
  wait "$pid" || { echo "ERROR: chunk processing failed (PID $pid)" >&2; exit 1; }
done
cat _chunk_*.out > "$NOT_PRESENT"
rm -f _chunk_* _chunk_*.out db_seq.txt tmp_query.tsv

end=$(date +%s)
runtime=$((end - start))
echo "Total runtime: ${runtime} seconds"
