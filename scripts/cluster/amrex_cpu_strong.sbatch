#!/usr/bin/env bash
# Strong MPI scaling; dependency on CPU validation is supplied by
# submit_campaign.sh, not hard-coded here. OpenMP gets a separate allocation
# so each thread count has the CPUs it actually claims.
#SBATCH --job-name=mhd-amrex-strong
set -euo pipefail
source "${MHD_CLUSTER_DIR:-$(cd "$(dirname "$0")" && pwd)}/lib_campaign.sh"
MHD_JOB_KIND=amrex-cpu-strong
campaign_init_environment
campaign_job_dir
campaign_build_dir
campaign_capture_environment
campaign_build_amrex_cpu
export OMP_PROC_BIND=close OMP_PLACES=cores
config="$(campaign_config_path "${MHD_CPU_CONFIG}")"
python3 "$MHD_SOURCE_DIR/scripts/scaling.py" \
    --executable "$MHD_BUILD_DIR/mhd2d" --config "$config" --mode mpi \
    --counts "$MHD_CPU_COUNTS" --repeats 5 --warmup 1 \
    --launcher 'srun --ntasks={p}' --output "$MHD_JOB_DIR/scaling_mpi.json" \
    2>&1 | tee "$MHD_JOB_DIR/scaling_mpi.log"
campaign_write_status pass "CPU strong MPI scaling completed"
