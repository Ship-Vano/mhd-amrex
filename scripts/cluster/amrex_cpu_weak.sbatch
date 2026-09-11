#!/usr/bin/env bash
# Weak MPI scaling.  The runner creates temporary scaled configs and preserves
# the summary, log and full environment in this job's immutable artifact.
#SBATCH --job-name=mhd-amrex-weak
set -euo pipefail
source "${MHD_CLUSTER_DIR:-$(cd "$(dirname "$0")" && pwd)}/lib_campaign.sh"
MHD_JOB_KIND=amrex-cpu-weak
campaign_init_environment
campaign_job_dir
campaign_build_dir
campaign_capture_environment
campaign_build_amrex_cpu
export OMP_NUM_THREADS=1 OMP_PROC_BIND=close OMP_PLACES=cores
config="$(campaign_config_path "${MHD_CPU_CONFIG}")"
python3 "$MHD_SOURCE_DIR/scripts/scaling.py" \
    --executable "$MHD_BUILD_DIR/mhd2d" --config "$config" --mode weak \
    --counts "$MHD_CPU_COUNTS" --weak-steps "$MHD_WEAK_STEPS" \
    --repeats 5 --warmup 1 --launcher 'srun --ntasks={p}' \
    --output "$MHD_JOB_DIR/scaling_weak.json" 2>&1 | tee "$MHD_JOB_DIR/scaling_weak.log"
campaign_write_status pass "CPU weak MPI scaling completed"
