#!/usr/bin/env bash
# One-GPU correctness gate, suitable for an RTX 4090 allocation.  This is not
# a performance benchmark: it first proves CUDA execution and CPU/GPU parity.
#SBATCH --job-name=mhd-amrex-cuda-validation
set -euo pipefail
source "${MHD_CLUSTER_DIR:-$(cd "$(dirname "$0")" && pwd)}/lib_campaign.sh"
MHD_JOB_KIND=amrex-cuda-validation
campaign_init_environment
campaign_job_dir
campaign_build_dir
campaign_capture_environment
campaign_run_cuda_validation
campaign_write_status pass "CUDA CTest and single-GPU CPU/GPU parity passed; no GPU performance claim"
