#!/usr/bin/env bash
# Corrected legacy smoke case.  The Python runner refuses non-clean/non-pinned
# input and performs the clone, patch, CTest and solver run under this artifact.
#SBATCH --job-name=mhd-legacy-corrected
set -euo pipefail
source "${MHD_CLUSTER_DIR:-$(cd "$(dirname "$0")" && pwd)}/lib_campaign.sh"
MHD_JOB_KIND=legacy-corrected
campaign_init_environment
campaign_require MHD_LEGACY_SOURCE
campaign_job_dir
campaign_capture_environment
python3 "$MHD_SOURCE_DIR/scripts/run_legacy_corrected.py" \
    --source "$MHD_LEGACY_SOURCE" --case "$MHD_LEGACY_CASE" \
    --mesh-backend structured --structured-nx "$MHD_LEGACY_NX" --structured-ny "$MHD_LEGACY_NY" \
    --artifact-dir "$MHD_JOB_DIR/legacy-run" --compiler "${CXX:-c++}" \
    --jobs "${MHD_BUILD_JOBS:-${SLURM_CPUS_PER_TASK:-8}}" \
    2>&1 | tee "$MHD_JOB_DIR/runner.log"
campaign_write_status pass "legacy_corrected isolated run completed"
