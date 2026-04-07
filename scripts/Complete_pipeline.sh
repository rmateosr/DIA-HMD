#!/bin/bash
# ABOUTME: Part of the DIA-NN Level 1 pipeline toolchain.
# ABOUTME: Submits DIA-NN search (Stage 1) then post-processing (Stage 2) with job dependency.
set -euo pipefail

SCRIPT_DIR="${SGE_O_WORKDIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")" && pwd)}}"
source "$SCRIPT_DIR/config.sh"

mkdir -p Library Reports log
chmod +x generate_diann_job.sh

./generate_diann_job.sh "$SAMPLE_DIR" "$FASTA_FILE" "$DIANN_IMG" > diann_search_job.sh
chmod +x diann_search_job.sh

DIANN_JOB=$(submit_job DIANN "$DIANN_THREADS" 6G "" diann_search_job.sh)
FILTER_JOB=$(submit_job StrictFilter 1 8G "$DIANN_JOB" strict_filter.sh)
# PostDIANN: 8 slots * 4G = 32G total — R analysis needs ~32 GB regardless of DIA-NN thread count
submit_job PostDIANN 8 4G "$FILTER_JOB" Post_DIANN_pipeline.sh
