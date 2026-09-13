# РП7 / T12 — CUDA-порт AMReX: состояние и первый runtime gate

Дата обновления: 2026-09-13.

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

Весь набор CTest проходит целиком (23/23) после переноса.

Главный вопрос к такому переносу — **изменил ли он ответ на CPU**. Проверено
прямо, а не понадеявшись на пороги тестов: собран решатель на ревизии до порта
(`74c1b3d`) теми же флагами, и оба бинарника прогнаны на всех девяти
конфигурациях из `inputs/`.

Совпали **точно**, до последней печатной цифры, на всех девяти случаях:

- число шагов и конечное время;
- счётчики `hlld_fallbacks`, `floor_events`, `nonpositive_cells` — важно,
  потому что порт переписал их с host-переменных на атомарные;
- диапазоны `rho` и `p`.

Разошлись две величины, и обе — на уровне округления:

| случай | `divB` max_abs до → после | нормированная |
|---|---|---|
| Брио–Ву | совпало побитово | совпало |
| постоянное состояние | совпало побитово | совпало |
| петля поля | 2.43e-16 → 2.54e-16 | 3.5e-15 → 3.7e-15 |
| Орзага–Танг | 4.48e-13 → 4.90e-13 | 5.2e-15 → 5.6e-15 |
| цилиндр | 1.38e-12 → 1.46e-12 | 4.8e-15 → 5.1e-15 |
| CP-альфвен | 1.46e-12 → 1.43e-12 | 2.9e-14 → 2.8e-14 |
| МГД-взрыв | 8.24e-13 → **4.19e-12** | 4.6e-15 → 2.3e-14 |

Пятикратный скачок у МГД-взрыва выглядит тревожно, но это свойство самой
величины, а не порта. `max|divB|` — максимум разности почти равных граневых
значений: оба операнда порядка единицы, их разность порядка `1e-16`, поэтому
изменение любого из них на один последний бит меняет разность на десятки
процентов. Такая величина не воспроизводима побитово между любыми двумя
сборками, отличающимися порядком операций, и её мантиссу нельзя цитировать как
результат. Содержательно заявляется другое: **нормированная норма остаётся на
уровне округления**, на два порядка ниже гейта `1e-12`, — и это выполняется во
всех случаях.

Дрейфы законов сохранения сдвинулись по той же причине (сумма накапливается в
другом порядке) и тоже остаются на уровне округления.

Вывод: перенос не изменил решение. Изменился только порядок операций, и
изменился он ровно там, где это ожидаемо и безвредно.

На рабочей станции Apple Silicon отсутствуют `nvcc` и NVIDIA driver, поэтому
CUDA compile/run здесь не выполнялся. Из этого не следует ни успешность, ни
неуспешность CUDA-порта: проверено, что порт не сломал CPU-путь, а не то, что
он работает на устройстве.

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
