#!/bin/bash
# ШАБЛОН: слабое масштабирование нового решателя. НЕ ПРОВЕРЕН НА ОЧЕРЕДИ.
# Заполнить <...> перед запуском; см. scripts/cluster/README.md.
#
# Слабое масштабирование отвечает на другой вопрос, чем сильное: не «во сколько
# раз быстрее», а «сколько стоит рост задачи». Число ячеек на ранг и число шагов
# фиксированы, домен растёт вместе с числом рангов, идеал -- постоянное время.
# Отклонение от постоянного времени -- это цена обменов и синхронизаций, и
# именно она интересна при переходе на несколько узлов.
#SBATCH --job-name=mhd2d-weak
#SBATCH --account=<ACCOUNT>
#SBATCH --partition=<PARTITION>
#SBATCH --nodes=<NODES>
#SBATCH --ntasks-per-node=<RANKS_PER_NODE>
#SBATCH --cpus-per-task=1
#SBATCH --time=<HH:MM:SS>
#SBATCH --output=slurm-%j-mhd2d-weak.out
set -euo pipefail

module purge
module load <COMPILER_MODULE> <MPI_MODULE> <CMAKE_MODULE>

ROOT="${SLURM_SUBMIT_DIR:-$PWD}"
BUILD="$ROOT/build/cluster-$SLURM_JOB_ID"
# Базовая сетка -- это сетка ОДНОГО ранга. При p рангах домен удваивается
# попеременно по x и по y, так что ячеек на ранг остаётся ровно столько же.
CONFIG="${CONFIG:-$ROOT/inputs/orszag_tang_uniform.json}"
OUT="$ROOT/benchmarks/raw/cluster/$SLURM_JOB_ID"
COUNTS="${COUNTS:-1,2,4,8,16,32}"
STEPS="${STEPS:-50}"
mkdir -p "$OUT"

cmake -S "$ROOT" -B "$BUILD" -DCMAKE_BUILD_TYPE=Release \
      -DMHD_MPI=ON -DMHD_OPENMP=OFF -DBUILD_TESTING=ON
cmake --build "$BUILD" --parallel 8
ctest --test-dir "$BUILD" --output-on-failure | tee "$OUT/ctest.log"

{
  echo "job_id=$SLURM_JOB_ID"
  echo "nodes=$SLURM_JOB_NUM_NODES"
  echo "ntasks=$SLURM_NTASKS"
  echo "nodelist=$SLURM_JOB_NODELIST"
  echo "commit=$(git -C "$ROOT" rev-parse HEAD)"
  echo "dirty=$(git -C "$ROOT" status --porcelain | wc -l)"
  echo "weak_steps=$STEPS"
  module list 2>&1
  lscpu 2>/dev/null || true
} > "$OUT/environment.txt"

export OMP_NUM_THREADS=1
export OMP_PROC_BIND=close
export OMP_PLACES=cores

python3 "$ROOT/scripts/scaling.py" \
    --executable "$BUILD/mhd2d" --config "$CONFIG" \
    --mode weak --counts "$COUNTS" --weak-steps "$STEPS" \
    --repeats 5 --warmup 1 \
    --launcher 'srun --ntasks={p}' \
    --output "$OUT/scaling_weak.json" | tee "$OUT/scaling_weak.log"

echo "результат: $OUT/scaling_weak.json"
