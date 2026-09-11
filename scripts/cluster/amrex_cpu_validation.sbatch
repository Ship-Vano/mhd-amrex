#!/usr/bin/env bash
# CPU correctness gate for one immutable campaign.  Submitted by
# submit_campaign.sh; variables come from its copied, user-owned site file.
#SBATCH --job-name=mhd-amrex-validate
set -euo pipefail
source "${MHD_CLUSTER_DIR:-$(cd "$(dirname "$0")" && pwd)}/lib_campaign.sh"
MHD_JOB_KIND=amrex-cpu-validation
campaign_init_environment
campaign_job_dir
campaign_build_dir
campaign_capture_environment
campaign_build_amrex_cpu
ctest --test-dir "$MHD_BUILD_DIR" --output-on-failure 2>&1 | tee "$MHD_JOB_DIR/ctest.log"
campaign_write_status pass "CPU Release build and CTest passed"
