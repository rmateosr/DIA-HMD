#!/bin/bash
# ABOUTME: Part of the DIA-NN Level 1 pipeline toolchain.
# ABOUTME: Generates diann_search_job.sh — emits a two-pass DIA-NN cluster job to stdout.
# Usage: ./generate_diann_job.sh /samples /fasta /diann.img > diann_search_job.sh

SAMPLE_DIR="$1"
FASTA_FILE="$2"
DIANN_IMG="$3"

cat <<'HEADER'
#!/bin/bash
set -o errexit
set -o nounset

SCRIPT_DIR="${SGE_O_WORKDIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")" && pwd)}}"
source "$SCRIPT_DIR/config.sh"
HEADER

# Pass 1: library-free search against Level 1 FASTA; generates predicted speclib
# Note: $DIANN_THREADS is escaped so it resolves at runtime (from config.sh), not at generation time.
# $FASTA_FILE is expanded now (baked into the generated script).
cat <<EOF

run_container /diann-2.0.2/diann-linux \\
--lib "" --threads \$DIANN_THREADS --verbose 1 \\
--out "Reports/report.parquet" \\
--qvalue 0.01 --matrices  --out-lib "Library/library.parquet" \\
--gen-spec-lib --predictor --fasta "$FASTA_FILE" \\
--fasta-search --min-fr-mz 200 --max-fr-mz 1800 --met-excision --min-pep-len 7 --max-pep-len 30 --min-pr-mz 300 --max-pr-mz 1800 --min-pr-charge 1 \\
--max-pr-charge 4 --cut K*,R* --missed-cleavages 1 --unimod4 --mass-acc 10 --mass-acc-ms1 4 --peptidoforms --reanalyse --rt-profiling --high-acc

EOF

# Pass 2: library-guided re-analysis using the predicted speclib from Pass 1
echo "run_container /diann-2.0.2/diann-linux \\"

for file in "$SAMPLE_DIR"/*.raw.dia; do
    echo "--f \"$file\" \\"
done

cat <<EOF
--lib "Library/library.predicted.speclib" \\
--threads \$DIANN_THREADS --verbose 1 --out "Reports/report_peptidoforms.tsv" \\
--qvalue 0.01 --matrices  --out-lib "Library/library_FROM_peptidoform.parquet" \\
--fasta "$FASTA_FILE" \\
--gen-spec-lib --met-excision --cut K*,R* --missed-cleavages 1 --unimod4 --mass-acc 10 --mass-acc-ms1 4.0 \\
--peptidoforms --reanalyse --rt-profiling --high-acc
EOF
