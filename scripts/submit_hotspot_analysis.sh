#!/bin/bash
# ABOUTME: Part of the DIA-NN Level 1 pipeline toolchain.
# ABOUTME: Runs hotspot peptide analysis in R.
[[ -n "${DEBUG:-}" ]] && set -xv
set -o errexit
set -o nounset
SCRIPT_DIR="${SGE_O_WORKDIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")" && pwd)}}"
source "$SCRIPT_DIR/config.sh"
ensure_tool Rscript R/4.4.3

# PROTEOME_FILE loaded from config.sh
Rscript noncanonicalpeptidesanalysis_Hotspot.R "$PROTEOME_FILE"


