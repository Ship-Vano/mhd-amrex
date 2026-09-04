#!/bin/bash
# ШАБЛОН: одна карта исторического решателя. НЕ ПРОВЕРЕН НА ОЧЕРЕДИ.
# Исторический решатель однопоточный: один ранг, одно ядро на задание.
#SBATCH --job-name=legacy-case
#SBATCH --account=<ACCOUNT>
#SBATCH --partition=<PARTITION>
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=<HH:MM:SS>
#SBATCH --output=slurm-%j-legacy-%x.out
set -euo pipefail

module purge
module load <COMPILER_MODULE> <CMAKE_MODULE> <PYTHON_MODULE>

ROOT="${SLURM_SUBMIT_DIR:-$PWD}"
SOURCE="${LEGACY_SOURCE:?set LEGACY_SOURCE to the immutable MHD2D checkout}"
CASE="${CASE:-rotor}"
NX="${NX:-128}"; NY="${NY:-128}"
OUT="$ROOT/benchmarks/raw/cluster/legacy-$CASE-$SLURM_JOB_ID"

# Раннер сам сверяет коммит источника и контрольную сумму архива и пишет
# манифест с quality gate даже при провале прогона.
python3 "$ROOT/scripts/run_legacy_corrected.py" \
    --source "$SOURCE" --case "$CASE" \
    --mesh-backend structured --structured-nx "$NX" --structured-ny "$NY" \
    --artifact-dir "$OUT" \
    --compiler "${CXX:-g++}" --jobs "${SLURM_CPUS_PER_TASK:-1}"

echo "manifest: $OUT/manifest.json"
