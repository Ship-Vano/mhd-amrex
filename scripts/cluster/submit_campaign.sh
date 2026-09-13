#!/usr/bin/env bash
# Submit an isolated T13 preparation campaign to SLURM.
# It only creates a new artifact directory and SLURM jobs; it never edits the
# mhd-amrex source checkout or the immutable legacy source.
set -euo pipefail

usage() {
    cat <<'EOF'
Usage:
  scripts/cluster/submit_campaign.sh --site scripts/cluster/sites/<site>.env \
      --legacy-source /path/to/MHD2D [--cuda-validation] [--dry-run]

The site file is copied into the new campaign directory.  Defaults are designed
for SLURM; use --dry-run first to inspect every sbatch command.
EOF
}

root="$(cd "$(dirname "$0")/../.." && pwd)"
site=""; legacy_override=""; cuda_validation=0; dry_run=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --site) site="${2:-}"; shift 2 ;;
        --legacy-source) legacy_override="${2:-}"; shift 2 ;;
        --cuda-validation|--gpu-probe) cuda_validation=1; shift ;;
        --dry-run) dry_run=1; shift ;;
        -h|--help) usage; exit 0 ;;
        *) echo "unknown argument: $1" >&2; usage >&2; exit 2 ;;
    esac
done
[[ -n "$site" && -f "$site" ]] || { echo "--site must name a readable site env file" >&2; exit 2; }
# shellcheck disable=SC1090
source "$site"
MHD_CLUSTER_DIR="$root/scripts/cluster"
MHD_SOURCE_DIR="${MHD_SOURCE_DIR:-$root}"
MHD_LEGACY_SOURCE="${legacy_override:-${MHD_LEGACY_SOURCE:-}}"
: "${MHD_ARTIFACT_ROOT:?MHD_ARTIFACT_ROOT is required in the site file}"
: "${MHD_LEGACY_SOURCE:?pass --legacy-source or set MHD_LEGACY_SOURCE in the site file}"
[[ -d "$MHD_SOURCE_DIR/.git" ]] || { echo "not a Git checkout: $MHD_SOURCE_DIR" >&2; exit 2; }
[[ -d "$MHD_LEGACY_SOURCE/.git" ]] || { echo "not a Git checkout: $MHD_LEGACY_SOURCE" >&2; exit 2; }
[[ -z "$(git -C "$MHD_SOURCE_DIR" status --porcelain=v1)" ]] || { echo "mhd-amrex source is dirty; commit or stash it before a reproducible campaign" >&2; exit 2; }
if (( ! dry_run )); then
    command -v sbatch >/dev/null || { echo "sbatch is unavailable; use run_ubuntu4090.sh for a non-SLURM host" >&2; exit 2; }
fi

for name in MHD_CPU_ACCOUNT MHD_CPU_PARTITION MHD_CPU_NODES MHD_CPU_RANKS_PER_NODE MHD_CPU_THREADS_PER_RANK MHD_CPU_TIME MHD_CPU_COUNTS MHD_OMP_COUNTS MHD_CPU_CONFIG MHD_WEAK_STEPS MHD_LEGACY_CASE MHD_LEGACY_NX MHD_LEGACY_NY; do
    [[ -n "${!name:-}" ]] || { echo "site file leaves $name unset" >&2; exit 2; }
done
max_count="$(tr ',' '\n' <<<"$MHD_CPU_COUNTS" | sort -n | tail -1)"
capacity=$((MHD_CPU_NODES * MHD_CPU_RANKS_PER_NODE))
(( max_count <= capacity )) || { echo "MHD_CPU_COUNTS asks for $max_count ranks, allocation has $capacity" >&2; exit 2; }
[[ "$MHD_CPU_THREADS_PER_RANK" = 1 ]] || { echo "this campaign's MPI strong/weak jobs are pure MPI; set MHD_CPU_THREADS_PER_RANK=1 and use MHD_OMP_COUNTS for OpenMP" >&2; exit 2; }
[[ -f "$( [[ "$MHD_CPU_CONFIG" = /* ]] && echo "$MHD_CPU_CONFIG" || echo "$MHD_SOURCE_DIR/$MHD_CPU_CONFIG" )" ]] || { echo "CPU config not found: $MHD_CPU_CONFIG" >&2; exit 2; }

MHD_EXPECTED_COMMIT="$(git -C "$MHD_SOURCE_DIR" rev-parse HEAD)"
commit="${MHD_EXPECTED_COMMIT:0:12}"
campaign="campaign-$(date -u +%Y%m%dT%H%M%SZ)-$commit"
campaign_dir="$MHD_ARTIFACT_ROOT/$campaign"
[[ ! -e "$campaign_dir" ]] || { echo "refusing to reuse campaign directory: $campaign_dir" >&2; exit 2; }
mkdir -p "$campaign_dir/jobs"
cp "$site" "$campaign_dir/site.env"
{
  echo "campaign_id=$campaign"
  echo "source_commit=$MHD_EXPECTED_COMMIT"
  echo "source_dirty_lines=0"
  echo "legacy_source=$MHD_LEGACY_SOURCE"
  echo "submitted_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "cuda_validation_requested=$cuda_validation"
} > "$campaign_dir/campaign.txt"

submit() {
    local label="$1"; shift
    local result
    if (( dry_run )); then
        printf 'DRY RUN %-18s' "$label" >&2; printf ' %q' "$@" >&2; printf '\n' >&2
        # Холостой прогон обещает показать команды, которые будут поданы, поэтому
        # он возвращает условный job id: иначе зависимости не попадут в вывод и
        # порядок заданий будет выглядеть иначе, чем при настоящей подаче.
        printf 'DRYRUN_JOBID_%s\n' "${label//-/_}"
        return
    fi
    result="$("$@")"
    printf '%-18s %s\n' "$label" "$result" | tee -a "$campaign_dir/submitted-jobs.txt" >&2
    # --parsable can append ;cluster. SLURM dependencies require only job ID.
    printf '%s\n' "${result%%;*}"
}

# Let SLURM inherit exported variables rather than serializing them in
# --export=NAME=value.  In particular, scaling lists contain commas, which are
# delimiters in SLURM's --export grammar and would otherwise be silently split.
export MHD_SOURCE_DIR MHD_LEGACY_SOURCE MHD_ARTIFACT_ROOT MHD_CAMPAIGN_ID MHD_EXPECTED_COMMIT MHD_CLUSTER_DIR
export MHD_MODULES MHD_CPU_CONFIG MHD_CPU_COUNTS MHD_OMP_COUNTS MHD_WEAK_STEPS
export MHD_LEGACY_CASE MHD_LEGACY_NX MHD_LEGACY_NY MHD_GPU_ARCH
exports="ALL"
common_cpu=(sbatch --parsable --export="$exports" --account="$MHD_CPU_ACCOUNT" --partition="$MHD_CPU_PARTITION" --nodes="$MHD_CPU_NODES" --ntasks-per-node="$MHD_CPU_RANKS_PER_NODE" --cpus-per-task="$MHD_CPU_THREADS_PER_RANK" --time="$MHD_CPU_TIME" --output="$campaign_dir/slurm-%j-%x.out")
# Замеры идут только после успешной валидации: считать быстро неверно смысла нет.
validation="$(submit amrex-validation "${common_cpu[@]}" "$root/scripts/cluster/amrex_cpu_validation.sbatch")"
cpu_after_validation=("${common_cpu[@]}" "--dependency=afterok:$validation")
submit amrex-strong "${cpu_after_validation[@]}" "$root/scripts/cluster/amrex_cpu_strong.sbatch"
submit amrex-weak "${cpu_after_validation[@]}" "$root/scripts/cluster/amrex_cpu_weak.sbatch"
max_omp="$(tr ',' '\n' <<<"$MHD_OMP_COUNTS" | sort -n | tail -1)"
omp=(sbatch --parsable --export="$exports" --account="$MHD_CPU_ACCOUNT" --partition="$MHD_CPU_PARTITION" --nodes=1 --ntasks=1 --cpus-per-task="$max_omp" --time="$MHD_CPU_TIME" --output="$campaign_dir/slurm-%j-%x.out")
omp+=("--dependency=afterok:$validation")
submit amrex-omp "${omp[@]}" "$root/scripts/cluster/amrex_cpu_omp.sbatch"
legacy_cpus="${MHD_LEGACY_BUILD_CPUS:-$MHD_CPU_RANKS_PER_NODE}"
legacy=(sbatch --parsable --export="$exports" --account="$MHD_CPU_ACCOUNT" --partition="$MHD_CPU_PARTITION" --nodes=1 --ntasks=1 --cpus-per-task="$legacy_cpus" --time="$MHD_CPU_TIME" --output="$campaign_dir/slurm-%j-%x.out")
submit legacy-corrected "${legacy[@]}" "$root/scripts/cluster/legacy_corrected_campaign.sbatch"
if (( cuda_validation )); then
    for name in MHD_GPU_ACCOUNT MHD_GPU_PARTITION MHD_GPU_NODES MHD_GPU_GPUS MHD_GPU_CPUS_PER_TASK MHD_GPU_TIME MHD_GPU_ARCH; do
        [[ -n "${!name:-}" ]] || { echo "--cuda-validation requires $name" >&2; exit 2; }
    done
    submit amrex-cuda-validation sbatch --parsable --export="$exports" --account="$MHD_GPU_ACCOUNT" --partition="$MHD_GPU_PARTITION" --nodes="$MHD_GPU_NODES" --ntasks=1 --cpus-per-task="$MHD_GPU_CPUS_PER_TASK" --gpus="$MHD_GPU_GPUS" --time="$MHD_GPU_TIME" --output="$campaign_dir/slurm-%j-%x.out" "$root/scripts/cluster/amrex_cuda_validation.sbatch" >/dev/null
fi
echo "campaign directory: $campaign_dir"
echo "download after completion: rsync -a <cluster>:${MHD_ARTIFACT_ROOT}/${campaign}/ ./$(basename "$campaign")/"
