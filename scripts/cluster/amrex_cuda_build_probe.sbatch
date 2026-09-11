#!/usr/bin/env bash
# Deprecated compatibility wrapper.  Use amrex_cuda_validation.sbatch; it runs
# the actual CUDA correctness gate and preserves parity logs.
#SBATCH --job-name=mhd-amrex-cuda-probe
set -euo pipefail
source "${MHD_CLUSTER_DIR:-$(cd "$(dirname "$0")" && pwd)}/lib_campaign.sh"
MHD_JOB_KIND=amrex-cuda-validation
campaign_init_environment
campaign_job_dir
campaign_build_dir
campaign_capture_environment
campaign_run_cuda_validation
campaign_write_status pass "CUDA CTest and single-GPU CPU/GPU parity passed; no GPU performance claim"
