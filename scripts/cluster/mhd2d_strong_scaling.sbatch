#!/bin/bash
# ШАБЛОН: сильное масштабирование нового решателя. НЕ ПРОВЕРЕН НА ОЧЕРЕДИ.
# Заполнить <...> перед запуском; см. scripts/cluster/README.md.
#
# Замеры делает scripts/scaling.py -- тот же код, что и на рабочей станции.
# Отдельной "кластерной" реализации замера нет намеренно: она разошлась бы с
# основной по определению медианы, числу повторов и проверке корректности, и
# кластерные числа стало бы не с чем сравнивать.
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
CONFIG="${CONFIG:-$ROOT/inputs/orszag_tang_uniform.json}"
OUT="$ROOT/benchmarks/raw/cluster/$SLURM_JOB_ID"
COUNTS="${COUNTS:-1,2,4,8,16,32}"
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

export OMP_NUM_THREADS=1
export OMP_PROC_BIND=close
export OMP_PLACES=cores

# 4. Сильное масштабирование по рангам MPI. Прогрев отбрасывается, повторов >= 5,
#    отчитывается медиана с MAD; прогон, изменивший диапазоны rho/p относительно
#    серийного, помечается как расхождение и скрипт возвращает ненулевой код.
python3 "$ROOT/scripts/scaling.py" \
    --executable "$BUILD/mhd2d" --config "$CONFIG" \
    --mode mpi --counts "$COUNTS" --repeats 5 --warmup 1 \
    --launcher 'srun --ntasks={p}' \
    --output "$OUT/scaling_mpi.json" | tee "$OUT/scaling_mpi.log"

# 5. То же по потокам OpenMP внутри одного узла (для сравнения с MPI).
if [ "${RUN_OMP:-1}" = "1" ]; then
  srun --ntasks=1 --cpus-per-task="${SLURM_CPUS_PER_TASK:-1}" \
    python3 "$ROOT/scripts/scaling.py" \
      --executable "$BUILD/mhd2d" --config "$CONFIG" \
      --mode omp --counts "${OMP_COUNTS:-1,2,4,8}" --repeats 5 --warmup 1 \
      --output "$OUT/scaling_omp.json" | tee "$OUT/scaling_omp.log"
fi

echo "результаты: $OUT/scaling_{mpi,omp}.json (схема совпадает с benchmarks/summary/)"
