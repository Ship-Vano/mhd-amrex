#!/usr/bin/env bash
# Isolated local campaign for Ubuntu, including a machine with RTX 4090.
# CPU validation and corrected legacy run are executed.  With --cuda-validation
# the script also performs the CUDA execution/parity gate on the allocated GPU.
set -euo pipefail

usage() { echo "Usage: $0 --legacy-source /path/to/MHD2D --artifact-root /path/to/artifacts [--cuda-validation]"; }
root="$(cd "$(dirname "$0")/../.." && pwd)"
legacy=""; artifact_root=""; cuda_validation=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --legacy-source) legacy="${2:-}"; shift 2 ;;
    --artifact-root) artifact_root="${2:-}"; shift 2 ;;
    --cuda-validation|--cuda-probe) cuda_validation=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
done
[[ -n "$legacy" && -d "$legacy/.git" ]] || { echo "--legacy-source must be a clean MHD2D Git checkout" >&2; exit 2; }
[[ -n "$artifact_root" ]] || { echo "--artifact-root is required" >&2; exit 2; }
command -v cmake >/dev/null; command -v python3 >/dev/null
[[ -z "$(git -C "$root" status --porcelain=v1)" ]] || { echo "mhd-amrex source is dirty; commit or stash it before a reproducible campaign" >&2; exit 2; }
MHD_EXPECTED_COMMIT="$(git -C "$root" rev-parse HEAD)"
commit="${MHD_EXPECTED_COMMIT:0:12}"
campaign="local-$(date -u +%Y%m%dT%H%M%SZ)-$commit"
MHD_SOURCE_DIR="$root" MHD_LEGACY_SOURCE="$legacy" MHD_ARTIFACT_ROOT="$artifact_root" MHD_CAMPAIGN_ID="$campaign" MHD_EXPECTED_COMMIT="$MHD_EXPECTED_COMMIT" MHD_MODULES=: MHD_JOB_KIND=amrex-cpu-validation \
  bash -c 'source "$MHD_SOURCE_DIR/scripts/cluster/lib_campaign.sh"; campaign_init_environment; campaign_job_dir; campaign_build_dir; campaign_capture_environment; campaign_build_amrex_cpu; ctest --test-dir "$MHD_BUILD_DIR" --output-on-failure 2>&1 | tee "$MHD_JOB_DIR/ctest.log"; campaign_write_status pass "CPU Release build and CTest passed"'
MHD_SOURCE_DIR="$root" MHD_LEGACY_SOURCE="$legacy" MHD_ARTIFACT_ROOT="$artifact_root" MHD_CAMPAIGN_ID="$campaign" MHD_EXPECTED_COMMIT="$MHD_EXPECTED_COMMIT" MHD_MODULES=: MHD_JOB_KIND=legacy-corrected \
  bash -c 'source "$MHD_SOURCE_DIR/scripts/cluster/lib_campaign.sh"; campaign_init_environment; campaign_job_dir; campaign_capture_environment; python3 "$MHD_SOURCE_DIR/scripts/run_legacy_corrected.py" --source "$MHD_LEGACY_SOURCE" --case rotor --mesh-backend structured --structured-nx 128 --structured-ny 128 --artifact-dir "$MHD_JOB_DIR/legacy-run" --compiler "${CXX:-c++}" --jobs "${MHD_BUILD_JOBS:-8}" 2>&1 | tee "$MHD_JOB_DIR/runner.log"; campaign_write_status pass "legacy_corrected isolated run completed"'
if (( cuda_validation )); then
  command -v nvidia-smi >/dev/null || { echo "--cuda-validation requested but no NVIDIA driver is visible" >&2; exit 2; }
  command -v nvcc >/dev/null || { echo "--cuda-validation requested but nvcc is unavailable" >&2; exit 2; }
  MHD_SOURCE_DIR="$root" MHD_ARTIFACT_ROOT="$artifact_root" MHD_CAMPAIGN_ID="$campaign" MHD_EXPECTED_COMMIT="$MHD_EXPECTED_COMMIT" MHD_MODULES=: MHD_GPU_ARCH=89 MHD_JOB_KIND=amrex-cuda-validation \
    bash -c 'source "$MHD_SOURCE_DIR/scripts/cluster/lib_campaign.sh"; campaign_init_environment; campaign_job_dir; campaign_build_dir; campaign_capture_environment; campaign_run_cuda_validation; campaign_write_status pass "CUDA CTest and single-GPU CPU/GPU parity passed; no GPU performance claim"'
fi
echo "campaign directory: $artifact_root/$campaign"
