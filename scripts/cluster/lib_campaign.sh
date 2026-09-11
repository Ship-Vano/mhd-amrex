#!/usr/bin/env bash
# Shared helpers for cluster jobs.  This file deliberately changes neither the
# source checkout nor an existing artifact directory.
set -euo pipefail

campaign_die() { echo "campaign: $*" >&2; exit 2; }

campaign_require() {
    local name="$1"
    [[ -n "${!name:-}" ]] || campaign_die "required variable $name is unset"
}

campaign_init_environment() {
    campaign_require MHD_SOURCE_DIR
    campaign_require MHD_ARTIFACT_ROOT
    campaign_require MHD_CAMPAIGN_ID
    [[ -d "$MHD_SOURCE_DIR/.git" ]] || campaign_die "MHD_SOURCE_DIR is not a Git checkout: $MHD_SOURCE_DIR"
    if [[ -n "${MHD_EXPECTED_COMMIT:-}" ]]; then
        local actual_commit
        actual_commit="$(git -C "$MHD_SOURCE_DIR" rev-parse HEAD)"
        [[ "$actual_commit" = "$MHD_EXPECTED_COMMIT" ]] ||
            campaign_die "source commit changed after submission ($actual_commit != $MHD_EXPECTED_COMMIT)"
        [[ -z "$(git -C "$MHD_SOURCE_DIR" status --porcelain=v1)" ]] ||
            campaign_die "source checkout became dirty after submission"
    fi
    # Site configuration is intentionally explicit.  It is sourced only from a
    # user-controlled file submitted with this campaign, never from job input.
    if [[ -n "${MHD_MODULES:-}" && "${MHD_MODULES}" != ":" ]]; then
        module purge
        eval "$MHD_MODULES"
    fi
    command -v cmake >/dev/null || campaign_die "cmake is unavailable after loading modules"
    command -v python3 >/dev/null || campaign_die "python3 is unavailable after loading modules"
}

campaign_job_dir() {
    campaign_require MHD_JOB_KIND
    local id="${SLURM_JOB_ID:-local-$(date -u +%Y%m%dT%H%M%SZ)}"
    MHD_JOB_DIR="$MHD_ARTIFACT_ROOT/$MHD_CAMPAIGN_ID/jobs/${MHD_JOB_KIND}-${id}"
    [[ ! -e "$MHD_JOB_DIR" ]] || campaign_die "refusing to overwrite job artifact: $MHD_JOB_DIR"
    mkdir -p "$MHD_JOB_DIR"
    export MHD_JOB_DIR
    # A scheduler records only an exit code.  Preserve an explicit status next
    # to logs even when `set -e` stops a job in the middle of a command.
    MHD_STATUS_WRITTEN=0
    trap campaign_record_failure EXIT
}

campaign_build_dir() {
    local job_id="${SLURM_JOB_ID:-local-$$}"
    local base="${MHD_SCRATCH_ROOT:-${SLURM_TMPDIR:-${TMPDIR:-/tmp}}}"
    MHD_BUILD_DIR="$base/mhd-amrex-${MHD_CAMPAIGN_ID}-${MHD_JOB_KIND}-${job_id}"
    [[ ! -e "$MHD_BUILD_DIR" ]] || campaign_die "refusing to reuse build directory: $MHD_BUILD_DIR"
    mkdir -p "$MHD_BUILD_DIR"
    export MHD_BUILD_DIR
}

campaign_capture_environment() {
    {
        echo "campaign_id=$MHD_CAMPAIGN_ID"
        echo "job_kind=$MHD_JOB_KIND"
        echo "job_id=${SLURM_JOB_ID:-local}"
        echo "submit_host=${SLURM_SUBMIT_HOST:-$(hostname)}"
        echo "nodes=${SLURM_JOB_NUM_NODES:-1}"
        echo "nodelist=${SLURM_JOB_NODELIST:-$(hostname)}"
        echo "ntasks=${SLURM_NTASKS:-1}"
        echo "cpus_per_task=${SLURM_CPUS_PER_TASK:-${OMP_NUM_THREADS:-1}}"
        echo "git_commit=$(git -C "$MHD_SOURCE_DIR" rev-parse HEAD)"
        echo "expected_commit=${MHD_EXPECTED_COMMIT:-not-enforced}"
        echo "git_dirty=$(git -C "$MHD_SOURCE_DIR" status --porcelain=v1 | wc -l | tr -d ' ')"
        echo "compiler=$(c++ --version 2>&1 | head -n 1 || true)"
        echo "cmake=$(cmake --version | head -n 1)"
        echo "python=$(python3 --version 2>&1)"
        echo "modules=$MHD_MODULES"
        command -v lscpu >/dev/null && lscpu || true
        command -v nvidia-smi >/dev/null && nvidia-smi -L || true
    } > "$MHD_JOB_DIR/environment.txt"
    git -C "$MHD_SOURCE_DIR" status --porcelain=v1 > "$MHD_JOB_DIR/source-status.txt"
    git -C "$MHD_SOURCE_DIR" rev-parse HEAD > "$MHD_JOB_DIR/source-commit.txt"
}

campaign_config_path() {
    local input="$1"
    if [[ "$input" = /* ]]; then echo "$input"; else echo "$MHD_SOURCE_DIR/$input"; fi
}

campaign_build_amrex_cpu() {
    cmake -S "$MHD_SOURCE_DIR" -B "$MHD_BUILD_DIR" \
        -DCMAKE_BUILD_TYPE=Release -DMHD_MPI=ON -DMHD_OPENMP=ON \
        -DMHD_HDF5=OFF -DBUILD_TESTING=ON 2>&1 | tee "$MHD_JOB_DIR/configure.log"
    cmake --build "$MHD_BUILD_DIR" --parallel "${MHD_BUILD_JOBS:-${SLURM_CPUS_PER_TASK:-8}}" 2>&1 | tee "$MHD_JOB_DIR/build.log"
}

campaign_build_amrex_cuda() {
    campaign_require MHD_GPU_ARCH
    cmake -S "$MHD_SOURCE_DIR" -B "$MHD_BUILD_DIR" \
        -DCMAKE_BUILD_TYPE=Release -DMHD_MPI=ON -DMHD_OPENMP=OFF \
        -DMHD_HDF5=OFF -DBUILD_TESTING=ON -DAMReX_GPU_BACKEND=CUDA \
        -DCMAKE_CUDA_ARCHITECTURES="$MHD_GPU_ARCH" -DAMReX_CUDA_FASTMATH=OFF \
        2>&1 | tee "$MHD_JOB_DIR/configure.log"
    cmake --build "$MHD_BUILD_DIR" --parallel "${MHD_BUILD_JOBS:-${SLURM_CPUS_PER_TASK:-8}}" 2>&1 | tee "$MHD_JOB_DIR/build.log"
}

campaign_build_amrex_cpu_reference() {
    cmake -S "$MHD_SOURCE_DIR" -B "$MHD_BUILD_DIR" \
        -DCMAKE_BUILD_TYPE=Release -DMHD_MPI=OFF -DMHD_OPENMP=OFF \
        -DMHD_HDF5=OFF -DBUILD_TESTING=OFF 2>&1 | tee "$MHD_JOB_DIR/cpu-reference-configure.log"
    cmake --build "$MHD_BUILD_DIR" --parallel "${MHD_BUILD_JOBS:-${SLURM_CPUS_PER_TASK:-8}}" 2>&1 | tee "$MHD_JOB_DIR/cpu-reference-build.log"
}

campaign_run_cuda_validation() {
    campaign_require MHD_GPU_ARCH
    command -v nvidia-smi >/dev/null || campaign_die "no NVIDIA driver/GPU visible in this allocation"
    command -v nvcc >/dev/null || campaign_die "nvcc is unavailable after loading GPU modules"
    nvidia-smi -L | tee "$MHD_JOB_DIR/nvidia-smi-L.txt"
    nvcc --version | tee "$MHD_JOB_DIR/nvcc-version.txt"

    local cuda_build="$MHD_BUILD_DIR"
    MHD_BUILD_DIR="$MHD_JOB_DIR/cpu-reference-build"
    campaign_build_amrex_cpu_reference
    MHD_BUILD_DIR="$cuda_build"
    campaign_build_amrex_cuda

    ctest --test-dir "$MHD_BUILD_DIR" -E '^mpi\\.decomposition_parity$' \
        --output-on-failure 2>&1 | tee "$MHD_JOB_DIR/ctest-cuda.log"
    python3 "$MHD_SOURCE_DIR/tests/check_cpu_gpu_parity.py" \
        --cpu "$MHD_JOB_DIR/cpu-reference-build/mhd2d" --gpu "$MHD_BUILD_DIR/mhd2d" \
        --config "$MHD_SOURCE_DIR/inputs/uniform_const.json" \
        --output-dir "$MHD_JOB_DIR/parity-uniform"
    python3 "$MHD_SOURCE_DIR/tests/check_cpu_gpu_parity.py" \
        --cpu "$MHD_JOB_DIR/cpu-reference-build/mhd2d" --gpu "$MHD_BUILD_DIR/mhd2d" \
        --config "$MHD_SOURCE_DIR/inputs/orszag_tang_uniform.json" \
        --output-dir "$MHD_JOB_DIR/parity-orszag-tang"
}

campaign_write_status() {
    local status="$1"; shift
    python3 - "$MHD_JOB_DIR/status.json" "$status" "$@" <<'PY'
import datetime, json, pathlib, sys
pathlib.Path(sys.argv[1]).write_text(json.dumps({
    "schema_version": 1, "status": sys.argv[2], "note": " ".join(sys.argv[3:]),
    "created_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
}, indent=2, sort_keys=True) + "\n")
PY
    MHD_STATUS_WRITTEN=1
}

campaign_record_failure() {
    local code=$?
    if (( code != 0 )) && [[ -n "${MHD_JOB_DIR:-}" ]] &&
       [[ "${MHD_STATUS_WRITTEN:-0}" != 1 ]]; then
        python3 - "$MHD_JOB_DIR/status.json" "$code" <<'PY' || true
import datetime, json, pathlib, sys
pathlib.Path(sys.argv[1]).write_text(json.dumps({
    "schema_version": 1, "status": "failed", "exit_code": int(sys.argv[2]),
    "created_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
}, indent=2, sort_keys=True) + "\n")
PY
    fi
}
