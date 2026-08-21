#!/bin/bash
# ABOUTME: CLI entry point for the DIANN paper pipeline.
# ABOUTME: Validates inputs, writes config, and runs the full pipeline.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

usage() {
    cat <<EOF
DIANN Paper Pipeline — Somatic mutation (hotspot) peptide detection from DIA-MS data

Usage: $(basename "$0") [OPTIONS]

Required:
  --input DIR         Path to directory containing *.raw.dia files
  --diann PATH        DIA-NN Apptainer image, Docker image name, or native binary path

Optional:
  --output DIR        Copy final results here (default: results/ in repo root)
  --fasta FILE        Custom FASTA file (default: bundled data/fasta/proteome.fasta)
  --proteome FILE     Non-mutated proteome to subtract (default: bundled data/fasta/proteome_nonmutated.fasta)
  --runtime RT        Container runtime for DIA-NN: apptainer, docker, native, auto (default: auto)
  --threads N         Number of threads for DIA-NN (default: 4)

Detection filters:
  --qvalue-gate Q     Per-run Q.Value a variant must reach to count as present (default: 0.001)
  --fragment-min F    Minimum share of library fragment intensity from ions containing the
                      mutated residue (default: 0.15)
  --min-replicates N  Injections of a sample a call must be seen in (default: 1 = off).
                      Needs --sample-map to know which runs belong together.

Cohort description (all optional):
  --truth FILE        Ground-truth TSV (Sample, Gene, Protein.Change, Detected.By.DIANN).
                      Without it the TP/FP classification stage is skipped.
  --sample-map FILE   TSV (run, sample, optional run_label) grouping injections into samples
  --aliases FILE      TSV (truth_name, run_name) for differing sample spellings
  --pools FILE        TSV (pool, members) naming samples that are mixtures of others
  --help              Show this help message

Examples:
  # Local machine with Docker
  $(basename "$0") --input /data/raw_files --diann biocontainers/diann:2.0.2 --runtime docker --threads 8

  # HPC with Apptainer image
  $(basename "$0") --input /data/raw_files --diann /apps/diann-2.0.2.img --threads 32

  # Native DIA-NN binary (no container)
  $(basename "$0") --input /data/raw_files --diann /usr/local/bin/diann-linux --runtime native

  # Reproduce the published cell-line cohort: replicate injections and a ground truth
  $(basename "$0") --input /data/celllines --diann diann-2.0.2.img --threads 32 \\
      --truth data/truth/Table1_hotspotcelllines_stopfree.tsv \\
      --sample-map data/cohorts/celllines_sample_map.tsv \\
      --aliases data/cohorts/celllines_aliases.tsv \\
      --pools data/cohorts/celllines_pools.tsv \\
      --min-replicates 2
EOF
    exit "${1:-0}"
}

# ---- Parse arguments ----
INPUT_DIR=""
DIANN_PATH=""
OUTPUT_DIR="$SCRIPT_DIR/results"
FASTA="$SCRIPT_DIR/data/fasta/proteome.fasta"
PROTEOME="$SCRIPT_DIR/data/fasta/proteome_nonmutated.fasta"
RUNTIME=""
THREADS=4
QVALUE_GATE=0.001
FRAGMENT_MIN=0.15
MIN_REPLICATES=1
TRUTH=""
SAMPLE_MAP=""
ALIASES=""
POOLS=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --input)    INPUT_DIR="$2"; shift 2 ;;
        --diann)    DIANN_PATH="$2"; shift 2 ;;
        --output)   OUTPUT_DIR="$2"; shift 2 ;;
        --fasta)    FASTA="$2"; shift 2 ;;
        --proteome) PROTEOME="$2"; shift 2 ;;
        --runtime)  RUNTIME="$2"; shift 2 ;;
        --threads)  THREADS="$2"; shift 2 ;;
        --qvalue-gate)    QVALUE_GATE="$2"; shift 2 ;;
        --fragment-min)   FRAGMENT_MIN="$2"; shift 2 ;;
        --min-replicates) MIN_REPLICATES="$2"; shift 2 ;;
        --truth)      TRUTH="$2"; shift 2 ;;
        --sample-map) SAMPLE_MAP="$2"; shift 2 ;;
        --aliases)    ALIASES="$2"; shift 2 ;;
        --pools)      POOLS="$2"; shift 2 ;;
        --help|-h)  usage 0 ;;
        *)          echo "ERROR: Unknown option '$1'" >&2; usage 1 ;;
    esac
done

# ---- Validate required arguments ----
if [[ -z "$INPUT_DIR" ]]; then
    echo "ERROR: --input is required" >&2
    usage 1
fi
if [[ -z "$DIANN_PATH" ]]; then
    echo "ERROR: --diann is required" >&2
    usage 1
fi
if [[ ! -d "$INPUT_DIR" ]]; then
    echo "ERROR: Input directory does not exist: $INPUT_DIR" >&2
    exit 1
fi
if [[ ! -f "$FASTA" ]]; then
    echo "ERROR: FASTA file not found: $FASTA" >&2
    exit 1
fi
if [[ ! -f "$PROTEOME" ]]; then
    echo "ERROR: Non-mutated proteome file not found: $PROTEOME" >&2
    exit 1
fi
# A named-but-missing cohort file is a typo, not "the generic case": fail rather than silently
# classifying against nothing or leaving replicates ungrouped.
for _pair in "truth:$TRUTH" "sample-map:$SAMPLE_MAP" "aliases:$ALIASES" "pools:$POOLS"; do
    _name="${_pair%%:*}"; _path="${_pair#*:}"
    if [[ -n "$_path" && ! -f "$_path" ]]; then
        echo "ERROR: --$_name file not found: $_path" >&2
        exit 1
    fi
done
if [[ "$MIN_REPLICATES" -gt 1 && -z "$SAMPLE_MAP" ]]; then
    echo "ERROR: --min-replicates $MIN_REPLICATES needs --sample-map to say which runs are" >&2
    echo "       injections of the same sample. Without it every run is its own sample and" >&2
    echo "       every call would be rejected." >&2
    exit 1
fi

# Check that input directory has .raw.dia files
shopt -s nullglob
RAW_FILES=("$INPUT_DIR"/*.raw.dia)
shopt -u nullglob
if [[ ${#RAW_FILES[@]} -eq 0 ]]; then
    echo "ERROR: No *.raw.dia files found in $INPUT_DIR" >&2
    exit 1
fi
echo "Found ${#RAW_FILES[@]} .raw.dia file(s) in $INPUT_DIR"

# ---- Write config ----
# Resolve to absolute paths
INPUT_DIR="$(cd "$INPUT_DIR" && pwd)"
FASTA="$(cd "$(dirname "$FASTA")" && pwd)/$(basename "$FASTA")"
PROTEOME="$(cd "$(dirname "$PROTEOME")" && pwd)/$(basename "$PROTEOME")"
for _var in TRUTH SAMPLE_MAP ALIASES POOLS; do
    _path="${!_var}"
    [[ -n "$_path" ]] && printf -v "$_var" '%s' "$(cd "$(dirname "$_path")" && pwd)/$(basename "$_path")"
done
# Resolve DIANN_PATH if it's a file (image or binary); leave Docker image names as-is
if [[ -e "$DIANN_PATH" ]]; then
    DIANN_PATH="$(cd "$(dirname "$DIANN_PATH")" && pwd)/$(basename "$DIANN_PATH")"
fi

# Resolve the container runtime and apptainer binary path now (on the login node,
# where the user has their environment set up) so scheduler jobs can use the
# absolute path directly — no dependency on module system layout or names.
APPTAINER_BIN=""
if [[ -z "$RUNTIME" || "$RUNTIME" == "auto" ]]; then
    if command -v apptainer &>/dev/null; then
        RUNTIME="apptainer"
        APPTAINER_BIN="$(command -v apptainer)"
    elif command -v singularity &>/dev/null; then
        RUNTIME="apptainer"
        APPTAINER_BIN="$(command -v singularity)"
    elif command -v docker &>/dev/null; then
        RUNTIME="docker"
    else
        RUNTIME="native"
    fi
else
    if [[ "$RUNTIME" == "apptainer" ]]; then
        APPTAINER_BIN="$(command -v apptainer 2>/dev/null || command -v singularity 2>/dev/null || true)"
    fi
fi

# Write runtime config as a separate override file (not sed-patching config.sh)
# so scheduler jobs can source it after run.sh exits.
RUN_CONFIG="$SCRIPT_DIR/scripts/config.run.sh"
cat > "$RUN_CONFIG" <<RUNEOF
# Auto-generated by run.sh -- do not edit; re-run run.sh to update.
SAMPLE_DIR="$INPUT_DIR"
FASTA_FILE="$FASTA"
DIANN_IMG="$DIANN_PATH"
PROTEOME_FILE="$PROTEOME"
CONTAINER_RUNTIME="$RUNTIME"
DIANN_THREADS=$THREADS
APPTAINER_BIN="${APPTAINER_BIN}"
QVALUE_GATE=$QVALUE_GATE
FRAGMENT_MIN=$FRAGMENT_MIN
MIN_REPLICATES=$MIN_REPLICATES
TRUTH_FILE="$TRUTH"
SAMPLE_MAP="$SAMPLE_MAP"
ALIASES_FILE="$ALIASES"
POOLS_FILE="$POOLS"
RUNEOF

echo ""
echo "=== Configuration ==="
echo "  Input:     $INPUT_DIR (${#RAW_FILES[@]} files)"
echo "  DIA-NN:    $DIANN_PATH"
echo "  FASTA:     $FASTA"
echo "  Proteome:  $PROTEOME"
echo "  Runtime:   ${RUNTIME:-auto}"
echo "  Threads:   $THREADS"
echo "  Filters:   Q.Value <= $QVALUE_GATE, fragment >= $FRAGMENT_MIN, $MIN_REPLICATES injection(s)"
if [[ -n "$TRUTH" ]]; then
    echo "  Truth:     $TRUTH"
else
    echo "  Truth:     none -- the TP/FP classification stage will be skipped"
fi
[[ -n "$SAMPLE_MAP" ]] && echo "  Samples:   $SAMPLE_MAP"
echo ""

# ---- Run pipeline ----
cd "$SCRIPT_DIR/scripts"
bash Complete_pipeline.sh

# ---- Copy results ----
if command -v sbatch &>/dev/null || command -v qsub &>/dev/null; then
    echo ""
    echo "=== Jobs submitted to scheduler ==="
    echo "Results will be in: $SCRIPT_DIR/scripts/Peptidomics_Results/"
    echo "After all jobs finish, copy results with:"
    echo "  cp -r $SCRIPT_DIR/scripts/Peptidomics_Results/* $OUTPUT_DIR/"
else
    if [[ -d "Peptidomics_Results" ]]; then
        mkdir -p "$OUTPUT_DIR"
        cp -r Peptidomics_Results/* "$OUTPUT_DIR/"
        echo ""
        echo "=== Results copied to $OUTPUT_DIR ==="
        ls -lh "$OUTPUT_DIR/"
    fi
fi
