#!/usr/bin/env bash
# Read-only preflight report for Ubuntu/RTX 4090 or a SLURM compute node.
# It never installs packages, changes the checkout, or probes the GPU with a
# numerical workload.  By default it prints to stdout; --output redirects all
# report output to an explicitly supplied path.
set -u

output=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --output) output="${2:-}"; shift 2 ;;
        -h|--help)
            echo "Usage: $0 [--output /path/to/preflight.txt]"
            exit 0
            ;;
        *) echo "unknown argument: $1" >&2; exit 2 ;;
    esac
done

if [[ -n "$output" ]]; then
    parent="$(dirname "$output")"
    mkdir -p "$parent"
    exec > "$output" 2>&1
fi

run() {
    local label="$1"; shift
    echo "--- $label ---"
    if command -v "$1" >/dev/null 2>&1; then
        "$@" || echo "[command failed: exit $?]"
    else
        echo "[not found: $1]"
    fi
    echo
}

echo "mhd-amrex preflight report"
echo "utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "host=$(hostname)"
echo "user=$(id -un)"
echo "cwd=$(pwd)"
echo

run "OS" bash -c 'if [[ -r /etc/os-release ]]; then cat /etc/os-release; else uname -a; fi'
run "kernel" uname -a
run "CPU topology" lscpu
run "CPU count" nproc
run "memory" free -h
run "filesystem" df -h .

echo "--- environment ---"
for name in CUDA_HOME CUDAToolkit_ROOT CUDA_PATH LD_LIBRARY_PATH PATH SLURM_JOB_ID SLURM_JOB_PARTITION SLURM_JOB_NODELIST SLURM_CPUS_PER_TASK SLURM_GPUS; do
    printf '%s=%s\n' "$name" "${!name-}"
done
echo

run "NVIDIA driver and GPU" nvidia-smi
run "NVIDIA GPU inventory" nvidia-smi --query-gpu=name,driver_version,memory.total,compute_cap --format=csv,noheader
run "CUDA compiler" nvcc --version
run "CMake" cmake --version
run "C++ compiler" c++ --version
run "C compiler" cc --version
run "Python" python3 --version
run "MPI launcher" mpirun --version
run "SLURM launcher" srun --version
run "SLURM submit" sbatch --version

echo "--- module environment ---"
if command -v module >/dev/null 2>&1; then
    module list 2>&1 || true
else
    echo "[module command not found]"
fi
echo

if [[ -d .git ]]; then
    echo "--- project checkout ---"
    git rev-parse --show-toplevel 2>/dev/null || true
    printf 'commit='; git rev-parse HEAD 2>/dev/null || true
    echo "status_porcelain:"
    git status --short 2>/dev/null || true
    echo
fi

echo "--- suggested interpretation ---"
echo "RTX 4090 should report compute capability 8.9; use CMAKE_CUDA_ARCHITECTURES=89."
echo "The first gate needs a working nvidia-smi, nvcc, CMake, C++ compiler and Python."
echo "For cluster jobs also provide scheduler/account/partition and MPI launcher details."
