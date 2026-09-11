# РП7 / T12 — CUDA-порт AMReX: состояние и первый runtime gate

Дата обновления: 2026-09-08.

## Что реализовано

`mhd2d` получил отдельный CUDA execution path на AMReX 25.01, не меняющий
математическую схему HLLD + MUSCL + CT + SSP-RK2.

- Все вычислительные обходы `MultiFab` в `src/MhdAmr.cpp` переведены с
  `LoopOnCpu` на `amrex::ParallelFor`/`amrex::For`; GPU не использует CPU
  tile-loop как вычислительный путь.
- Постановка задачи теперь `DeviceProblem`: trivially-copyable POD, захватываемый
  device-lambda по значению. `std::function` и захват `this` из GPU-пути
  исключены.
- Физические `ext_dir`-границы работают через `GpuBndryFuncFab<ExtDirGpuFill>`.
- Счётчики HLLD→HLL и срабатывания floor собираются атомарно в
  `Gpu::DeviceScalar`; `ComputeDt`, `MaxDivB`, энергия и диапазоны используют
  `ReduceOps` вместо ссылок на host accumulator.
- Gas-flux reflux при AMR выбирает `RunOn::Gpu` в GPU-сборке. CT и порядок
  стадий сохранены.
- CMake вызывает `setup_target_for_cuda_compilation(mhd2d)` только для
  `AMReX_GPU_BACKEND=CUDA`; preset `cuda-release` выставляет
  `CMAKE_CUDA_ARCHITECTURES=89` для RTX 4090 и отключает CUDA fast math для
  осмысленного первого CPU/GPU parity gate.

Архитектурный тест `arch.gpu_portability` выполняется даже на CPU-host: он
ловит возвращение `LoopOnCpu`, `std::function` в problem data, потерю GPU
boundary functor, device counters/reductions или CUDA setup. Это статический
gate, не доказательство запуска на GPU.

## Что уже проверено локально

CPU Release пересобран после переноса. Пройдены:

```text
solver.compute_only
canonical.constant_state
arch.gpu_portability
```

На рабочей станции Apple Silicon отсутствуют `nvcc` и NVIDIA driver, поэтому
CUDA compile/run здесь не выполнялся. Из этого не следует ни успешность, ни
неуспешность CUDA-порта.

## Первый gate на RTX 4090

Нужны CUDA Toolkit 12.x, совместимый NVIDIA driver, CMake, GCC/G++ и Python.
Сначала выполняется именно correctness gate, затем — performance work:

```sh
cmake --preset cpu-release && cmake --build --preset cpu-release
cmake --preset cuda-release && cmake --build --preset cuda-release
ctest --preset cuda-release -E '^mpi\.decomposition_parity$'
python3 tests/check_cpu_gpu_parity.py \
  --cpu build/cpu-release/mhd2d --gpu build/cuda-release/mhd2d \
  --config inputs/uniform_const.json --output-dir benchmarks/raw/cuda/parity-uniform
python3 tests/check_cpu_gpu_parity.py \
  --cpu build/cpu-release/mhd2d --gpu build/cuda-release/mhd2d \
  --config inputs/orszag_tang_uniform.json --output-dir benchmarks/raw/cuda/parity-orszag
```

Parity требует совпадения счётчиков (`fallbacks`, `floors`, `nonpositive`) и
сопоставляет `rho/p` и нормы `div B` с `atol=5e-11`, `rtol=5e-10`. Сырые логи
двух бинарников записываются в `--output-dir`, чтобы их можно было передать
анализатору без повторного запуска.

Для Ubuntu есть одна команда, создающая изолированный campaign artifact:

```sh
scripts/cluster/run_ubuntu4090.sh --legacy-source /path/to/MHD2D \
  --artifact-root /data/mhd-artifacts --cuda-validation
```

Для SLURM/K10 добавляется `--cuda-validation` к
`scripts/cluster/submit_campaign.sh`. Job записывает версии `nvcc`, GPU,
configure/build/CTest и parity logs. Сборки создаются в artifact/scratch;
source checkout остаётся только для чтения.

## Чего этот change set не утверждает

- Нет выполненного CUDA runtime/parity результата до запуска на RTX 4090.
- Нет multi-rank MPI + GPU gate, GPU scaling, profiling Nsight или сравнения
  производительности с CPU.
- Нет основания заявлять ускорение, полную GPU-портируемость ввода-вывода или
  завершённость T12. Это следующие gate после успешного single-GPU parity.
