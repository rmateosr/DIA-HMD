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

# Lay the proteome out one protein per line so peptide sequences can be found with grep -F.
# One line per protein and not one string for the whole file: a grep -F match may not run past
# the end of a line, so a peptide is only ever compared against a sequence that exists. Joining
# the proteins would let a peptide straddling the boundary between two of them match a chimera
# no protein produces, and be discarded as canonical on that evidence.
# Records are cut on lines beginning with ">" rather than on the character: twelve UniProt
# descriptions contain one ("DNA dC->dU-editing enzyme"), and splitting there would break those
# proteins in two, which is the same error in the opposite direction.
# Ile and Leu have identical elemental composition and are unresolvable by MS, so both
# operands are collapsed to a single letter: the comparison asks whether a peptide is
# indistinguishable by mass from something canonical, not whether it matches literally.
awk '/^>/ {if (NR>1) printf "\n"; next} {printf "%s", $0} END {printf "\n"}' "$DB" | tr 'I' 'L' > db_seq.txt

# Convert FASTA to header<TAB>sequence TSV
awk 'BEGIN{RS=">"; ORS=""} NR>1 {n=split($0, lines, "\n"); header=lines[1]; seq=""; for (i=2; i<=n; i++) seq=seq lines[i]; print header "\t" seq "\n"}' "$QUERY" > tmp_query.tsv

> "$NOT_PRESENT"

# Split input into chunks and search in parallel using background jobs
NJOBS=8
split -n l/$NJOBS tmp_query.tsv _chunk_
pids=()
for chunk in _chunk_*; do
  (while IFS=$'\t' read -r header seq; do
    if ! grep -m 1 -qF "$(printf '%s' "$seq" | tr 'I' 'L')" db_seq.txt; then
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
