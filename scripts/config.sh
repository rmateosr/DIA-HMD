#!/bin/bash
# ABOUTME: Central configuration for the DIANN paper pipeline.
# ABOUTME: Edit the values below before running Complete_pipeline.sh.

# ---- User-configurable paths ----
SAMPLE_DIR="/path/to/your/DIA/raw/files"           # directory containing *.raw.dia
FASTA_FILE="../data/fasta/proteome.fasta"           # bundled with this repo
DIANN_IMG="/path/to/diann-2.0.2.img"                # DIA-NN Apptainer/Docker image or native binary path
PROTEOME_FILE="../data/fasta/proteome_nonmutated.fasta"  # bundled with this repo

# ---- Detection filters ----
# Per-run Q.Value a variant precursor must reach to count as present in that sample. DIA-NN's
# own output is already at 1% run FDR, so this is the gate that actually decides presence.
QVALUE_GATE=0.001
# Share of a variant peptide's library fragment intensity that must come from ions containing the
# mutated residue. A fragment that misses the site has the same mass in the mutant as in the wild
# type, so it is signal the normal protein produces identically.
FRAGMENT_MIN=0.15
# Injections of a sample a call must be seen in. 1 = off, which is the default because your
# replicate structure is unknown to this pipeline. The published cohorts used 2, which needs a
# SAMPLE_MAP below to say which runs belong together.
MIN_REPLICATES=1

# ---- Cohort description (all optional; leave empty for the generic case) ----
# Ground truth for TP/FP classification: TSV with Sample, Gene, Protein.Change,
# Detected.By.DIANN. Empty means the classification stage is skipped -- there is nothing to
# classify against. See data/truth/ for the two published cohorts.
TRUTH_FILE=""
# TSV (run, sample, optional run_label) grouping injections into samples. Empty: one sample per
# run. See data/cohorts/.
SAMPLE_MAP=""
# TSV (truth_name, run_name) for samples the truth table spells differently from the run names.
ALIASES_FILE=""
# TSV (pool, members) naming samples that are mixtures of others.
POOLS_FILE=""

# ---- Container runtime for DIA-NN ----
# Options: "apptainer", "docker", "native", or "" for auto-detect
# - apptainer: apptainer exec $DIANN_IMG /diann-2.0.2/diann-linux ...
# - docker:    docker run --rm -v ... $DIANN_IMG /diann-2.0.2/diann-linux ...
# - native:    runs DIA-NN binary directly (DIANN_IMG must be path to the binary)
CONTAINER_RUNTIME=""

# ---- DIA-NN threads (used by generate_diann_job.sh) ----
DIANN_THREADS=8  # increase if you have more cores available

# ---- Module system (leave empty if tools are already on PATH) ----
MODULE_BASE=""  # e.g., "/usr/local/package/modulefiles/"
APPTAINER_BIN=""  # resolved at run.sh time; overridden by config.run.sh

# Override defaults with values written by run.sh (keeps this template untouched)
_config_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ -f "${_config_dir}/config.run.sh" ]]; then
  source "${_config_dir}/config.run.sh"
fi

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

# ---- Helper: stop with a readable message when an upstream stage produced nothing ----
# SGE's -hold_jid releases a job when its predecessor FINISHES, whatever its exit status, unlike
# SLURM's --dependency=afterok. So on SGE a failed DIA-NN search does not stop the chain: the next
# stage starts and dies on a missing file, and what the user sees is a stack trace from whichever
# tool happened to open it first. These checks turn that into one line naming the stage to look at.
# Usage: require_inputs <upstream job name> <file>...
require_inputs() {
  local upstream="$1"; shift
  local missing=0 f
  for f in "$@"; do
    if [[ ! -s "$f" ]]; then
      echo "ERROR: required input is missing or empty: $f" >&2
      missing=1
    fi
  done
  if (( missing )); then
    echo "ERROR: the $upstream stage did not produce what this stage needs." >&2
    echo "       Look at log/${upstream}.* for why, then rerun from there." >&2
    exit 1
  fi
}

# ---- Helper: directories the container has to be able to see ----
# The ones we name, plus wherever any symlink among the inputs actually points. Keeping raw data
# elsewhere and symlinking it into an input directory is a normal way to organise it, and binding
# the directory does not bind what its links point at.
#
# Whether that matters depends on the runtime. Docker mounts nothing but what -v names, so a
# symlink out of the input directory always dangles inside the container. Apptainer is
# configuration-dependent: a site with a permissive bind path list may expose the target anyway,
# and one with a restrictive list will not. Resolving the inputs and binding their real parents
# makes the outcome the same everywhere, and costs nothing when there are no symlinks -- readlink
# -f returns the path itself and the duplicate is dropped.
_container_bind_dirs() {
  local dirs=("$PWD" "$SAMPLE_DIR" "$(dirname "$FASTA_FILE")")
  local f real
  for f in "$SAMPLE_DIR"/*.raw.dia "$FASTA_FILE"; do
    [[ -e "$f" ]] || continue
    real="$(readlink -f "$f" 2>/dev/null)" || continue
    [[ -n "$real" ]] && dirs+=("$(dirname "$real")")
  done
  printf '%s\n' "${dirs[@]}" | awk 'NF && !seen[$0]++'
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
  local apptainer_cmd="${APPTAINER_BIN:-}"
  if [[ -z "$apptainer_cmd" ]]; then
    apptainer_cmd="$(command -v apptainer 2>/dev/null || command -v singularity 2>/dev/null || true)"
  fi
  if [[ "$runtime" == "apptainer" && -z "$apptainer_cmd" ]]; then
    echo "ERROR: apptainer/singularity not found. Run run.sh from a node where apptainer is on PATH." >&2
    exit 1
  fi
  case "$runtime" in
    apptainer)
      local bind_args=() d
      while IFS= read -r d; do bind_args+=(--bind "$d"); done < <(_container_bind_dirs)
      "$apptainer_cmd" exec "${bind_args[@]}" "$DIANN_IMG" "$@"
      ;;
    docker)
      local vol_args=() d
      while IFS= read -r d; do vol_args+=(-v "$d:$d"); done < <(_container_bind_dirs)
      docker run --rm "${vol_args[@]}" -w "$PWD" "$DIANN_IMG" "$@"
      ;;
    native)
      # $1 is the container-internal command path — not meaningful in native mode
      shift
      local diann_dir
      diann_dir="$(cd "$(dirname "$DIANN_IMG")" && pwd)"
      # Create libgomp symlink if missing (DIA-NN bundles it with a hashed name)
      if [[ ! -e "${diann_dir}/libgomp.so.1" ]] && ls "${diann_dir}"/libgomp-*.so.1 &>/dev/null; then
          ln -sf "${diann_dir}"/libgomp-*.so.1 "${diann_dir}/libgomp.so.1"
      fi
      export LD_LIBRARY_PATH="${diann_dir}:${LD_LIBRARY_PATH:-}"
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
    local rc=${PIPESTATUS[0]}
    if [[ $rc -ne 0 ]]; then
        echo "ERROR: $name failed with exit code $rc. See log/${name}.log" >&2
        return $rc
    fi
    echo "local_${name}"
  fi
}
