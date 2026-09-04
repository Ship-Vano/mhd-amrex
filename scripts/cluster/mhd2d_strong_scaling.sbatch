#!/bin/bash
# ШАБЛОН: сильное масштабирование нового решателя. НЕ ПРОВЕРЕН НА ОЧЕРЕДИ.
# Заполнить <...> перед запуском; см. scripts/cluster/README.md.
#SBATCH --job-name=mhd2d-strong
#SBATCH --account=<ACCOUNT>
#SBATCH --partition=<PARTITION>
#SBATCH --nodes=<NODES>
#SBATCH --ntasks-per-node=<RANKS_PER_NODE>
#SBATCH --cpus-per-task=<THREADS_PER_RANK>
#SBATCH --time=<HH:MM:SS>
#SBATCH --output=slurm-%j-mhd2d-strong.out
set -euo pipefail

module purge
module load <COMPILER_MODULE> <MPI_MODULE> <CMAKE_MODULE>

ROOT="${SLURM_SUBMIT_DIR:-$PWD}"
BUILD="$ROOT/build/cluster-$SLURM_JOB_ID"
CONFIG="${CONFIG:-$ROOT/inputs/orszag_tang.json}"
OUT="$ROOT/benchmarks/raw/cluster/$SLURM_JOB_ID"
mkdir -p "$OUT"

# 1. Собираем там же, где считаем.
cmake -S "$ROOT" -B "$BUILD" -DCMAKE_BUILD_TYPE=Release \
      -DMHD_MPI=ON -DMHD_OPENMP=ON -DBUILD_TESTING=ON
cmake --build "$BUILD" --parallel "${SLURM_CPUS_PER_TASK:-8}"

# 2. Регрессии до замеров: считать быстро неверно смысла нет.
ctest --test-dir "$BUILD" --output-on-failure | tee "$OUT/ctest.log"

# 3. Окружение прогона -- в манифест.
{
  echo "job_id=$SLURM_JOB_ID"
  echo "nodes=$SLURM_JOB_NUM_NODES"
  echo "ntasks=$SLURM_NTASKS"
  echo "cpus_per_task=${SLURM_CPUS_PER_TASK:-1}"
  echo "nodelist=$SLURM_JOB_NODELIST"
  echo "commit=$(git -C "$ROOT" rev-parse HEAD)"
  echo "dirty=$(git -C "$ROOT" status --porcelain | wc -l)"
  module list 2>&1
  lscpu 2>/dev/null || true
} > "$OUT/environment.txt"

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export OMP_PROC_BIND=close
export OMP_PLACES=cores

# 4. Развёртка по числу рангов. Прогрев отбрасывается, повторов >= 5.
for RANKS in 1 2 4 8 16 32; do
    [ "$RANKS" -le "$SLURM_NTASKS" ] || continue
    srun --ntasks="$RANKS" "$BUILD/mhd2d" "$CONFIG" > "$OUT/warmup_$RANKS.log" 2>&1 || true
    for REP in 1 2 3 4 5; do
        /usr/bin/time -f "%e" -o "$OUT/time_${RANKS}_${REP}.txt" \
            srun --ntasks="$RANKS" "$BUILD/mhd2d" "$CONFIG" \
            > "$OUT/run_${RANKS}_${REP}.log" 2>&1
    done
done

echo "raw results in $OUT; summarise with scripts/scaling.py conventions"
