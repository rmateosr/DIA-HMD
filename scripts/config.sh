#!/bin/bash
# ABOUTME: Central configuration for the DIANN paper pipeline.
# ABOUTME: Edit the values below before running Complete_pipeline.sh.

# ---- User-configurable paths ----
SAMPLE_DIR="/path/to/your/DIA/raw/files"           # directory containing *.raw.dia
FASTA_FILE="../data/fasta/proteome.fasta"           # bundled with this repo
DIANN_IMG="/path/to/diann-2.0.2.img"                # DIA-NN Apptainer/Docker image or native binary path
PROTEOME_FILE="../data/fasta/human_canonical_proteome.fasta"  # bundled with this repo

# ---- Container runtime for DIA-NN ----
# Options: "apptainer", "docker", "native", or "" for auto-detect
# - apptainer: apptainer exec $DIANN_IMG /diann-2.0.2/diann-linux ...
# - docker:    docker run --rm -v ... $DIANN_IMG /diann-2.0.2/diann-linux ...
# - native:    runs DIA-NN binary directly (DIANN_IMG must be path to the binary)
CONTAINER_RUNTIME=""

# ---- DIA-NN threads (used by generate_diann_job.sh) ----
DIANN_THREADS=4  # increase if you have more cores available

# ---- Module system (leave empty if tools are already on PATH) ----
MODULE_BASE=""  # e.g., "/usr/local/package/modulefiles/"

# ---- Helper: load a tool via module system only if not already on PATH ----
ensure_tool() {
  local tool="$1" module_name="$2"
  if command -v "$tool" &>/dev/null; then return 0; fi
  if type module &>/dev/null && [[ -n "${MODULE_BASE:-}" ]]; then
    module use "$MODULE_BASE"
    module load "$module_name"
  else
    echo "ERROR: '$tool' not found on PATH and no module system available." >&2
    echo "Install $tool or set MODULE_BASE in config.sh" >&2
    exit 1
  fi
}

# ---- Helper: run a command inside the DIA-NN container (or natively) ----
# Usage: run_container <command> [args...]
# The container image is taken from $DIANN_IMG. The runtime is auto-detected
# unless CONTAINER_RUNTIME is set.
run_container() {
  local runtime="${CONTAINER_RUNTIME:-auto}"
  if [[ "$runtime" == "auto" ]]; then
    if command -v apptainer &>/dev/null; then runtime="apptainer"
    elif command -v singularity &>/dev/null; then runtime="apptainer"
    elif command -v docker &>/dev/null; then runtime="docker"
    else runtime="native"; fi
  fi
  case "$runtime" in
    apptainer)
      apptainer exec "$DIANN_IMG" "$@"
      ;;
    docker)
      docker run --rm \
        -v "$PWD:$PWD" \
        -v "$(dirname "$FASTA_FILE"):$(dirname "$FASTA_FILE")" \
        -v "$SAMPLE_DIR:$SAMPLE_DIR" \
        -w "$PWD" \
        "$DIANN_IMG" "$@"
      ;;
    native)
      # In native mode, DIANN_IMG should be the path to the diann-linux binary
      "$DIANN_IMG" "$@"
      ;;
    *)
      echo "ERROR: Unknown CONTAINER_RUNTIME '$runtime'. Use apptainer, docker, or native." >&2
      exit 1
      ;;
  esac
}

# ---- Helper: submit a job to SGE, SLURM, or run locally ----
# Usage: JOB_ID=$(submit_job <name> <slots> <mem_per_slot> <hold_jid|""> <script>)
# Returns the job ID on stdout (capture it for dependency chains).
submit_job() {
  local name="$1" slots="$2" mem_per_slot="$3" hold_jid="$4" script="$5"

  if command -v sbatch &>/dev/null; then
    # SLURM
    local total_mem=$(( ${mem_per_slot%G} * slots ))
    local args="--job-name=$name --cpus-per-task=$slots --mem=${total_mem}G"
    args+=" --output=log/%x_%j.out --error=log/%x_%j.err"
    [[ -n "$hold_jid" ]] && args+=" --dependency=afterok:$hold_jid"
    sbatch $args --parsable "$script"
  elif command -v qsub &>/dev/null; then
    # SGE
    local args="-N $name -pe def_slot $slots -l s_vmem=${mem_per_slot} -cwd -o log -e log -S /bin/bash"
    [[ -n "$hold_jid" ]] && args+=" -hold_jid $hold_jid"
    qsub $args "$script" | grep -oP 'Your job \K\d+'
  else
    # Local fallback — run sequentially
    echo "[local] Running $name..." >&2
    bash "$script" 2>&1 | tee "log/${name}.log"
    echo "local_${name}"
  fi
}
