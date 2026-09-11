#!/usr/bin/env bash
# Strong OpenMP scaling on exactly one allocated rank. Kept separate from MPI
# scaling so SLURM CPU affinity matches OMP_NUM_THREADS for every repeat.
#SBATCH --job-name=mhd-amrex-omp
set -euo pipefail
source "${MHD_CLUSTER_DIR:-$(cd "$(dirname "$0")" && pwd)}/lib_campaign.sh"
MHD_JOB_KIND=amrex-cpu-omp
campaign_init_environment
campaign_job_dir
campaign_build_dir
campaign_capture_environment
campaign_build_amrex_cpu
config="$(campaign_config_path "${MHD_CPU_CONFIG}")"
python3 "$MHD_SOURCE_DIR/scripts/scaling.py" \
    --executable "$MHD_BUILD_DIR/mhd2d" --config "$config" --mode omp \
    --counts "$MHD_OMP_COUNTS" --repeats 5 --warmup 1 \
    --output "$MHD_JOB_DIR/scaling_omp.json" 2>&1 | tee "$MHD_JOB_DIR/scaling_omp.log"
campaign_write_status pass "CPU strong OpenMP scaling completed"
