# Как запускать: сборка, тесты, сценарии, кластер

Один документ на оба решателя. Всё, что здесь описано, воспроизводится с
чистого клона; команды даны целиком, без подстановок «по смыслу».

---

## 0. Что где лежит

| | новый решатель | исторический решатель |
|---|---|---|
| профиль | `mhd2d_amrex` | `legacy_corrected` |
| исходники | `src/` (этот репозиторий) | `/Users/ivansamanov/Documents/MHD2D`, **не изменяется** |
| как получается | собирается напрямую | клон @ `9d0f60e` + `legacy/patches/0001-legacy-corrected-physics.patch` |
| сетка | декартова, AMR через AMReX | неструктурированная треугольная |
| запуск | `./build/release/mhd2d <config.json>` | `scripts/run_legacy_corrected.py` |

Исторический снимок закреплён тегом `legacy_vkr/9d0f60e` и контрольной суммой
архива (`docs/LEGACY_BASELINE.md`); раннер сверяет их при каждом запуске и
отказывается работать при расхождении или грязном дереве.

---

## 1. Сборка

```sh
cmake --preset release          # MPI + OpenMP, Release
cmake --build --preset release -j 8
ctest --preset release          # 22 теста, ~15 c
```

Другие пресеты (`cmake --list-presets`):

| пресет | зачем |
|---|---|
| `release` | обычная работа: MPI + OpenMP |
| `cpu-release` | серийный, для чистых замеров |
| `cpu-debug` | `AMReX_ASSERTIONS` + `AMReX_BOUND_CHECK` — гонять перед коммитом ядра |
| `mpi-release` | MPI без OpenMP: масштабирование по рангам |
| `profile` | `AMReX_TINY_PROFILE` — профиль горячего пути |
| `hdf5-release` | вывод HDF5 |
| `cuda-release` | CUDA Release для NVIDIA; preset нацелен на Ada `sm_89` (RTX 4090) |

---

## 2. Тесты: что именно проверяется

```sh
ctest --preset release                 # всё
ctest --preset release -R canonical    # только канонические задачи
ctest --preset release -R amr          # только консервативность AMR
ctest --preset release --output-on-failure -R briowu.independent_reference
```

| тест | что упадёт, если сломать |
|---|---|
| `kernel.unit`, `kernel.numerics` | инварианты HLLD, лимитеров, SSP-RK2, одномерный предел угловой ЭДС |
| `standalone.briowu`, `.alfven32`, `.alfven_order` | регрессии схемы; наблюдаемый порядок $\geq1.8$ |
| `standalone.dai_woodward`, `.loop` | вторая римановская задача и перенос петли поля |
| `canonical.constant_state` | постоянное состояние обязано сохраняться **побитово** |
| `canonical.orszag_tang`, `.rotor` | диапазоны `rho`/`p` против литературных цветовых шкал |
| `canonical.mhd_blast` | низкое `beta`: откат HLLD→HLL срабатывает, состояние остаётся положительным |
| `amr.conservation` | согласование потоков на стыке уровней; тест сам падает, если стыка нет или если при `reflux=false` дефект не проявился |
| `amr.regrid_conservation` | перестроение сетки не двигает интегралы |
| `mpi.decomposition_parity` | результат не зависит от числа рангов (1/2/4) |
| `briowu.independent_reference` | согласие с независимой схемой (Куртганова–Тадмора) |
| `arch.kernel_purity` | слой ядер не тянет контейнеры AMReX (ADR 0001) |
| `arch.gpu_portability` | статически ловит host-only loop/capture и потерю CUDA boundary/reduction пути |
| `config.*` | строгая схема конфигурации, режимы аблации реально влияют |
| `manifest.repeatable` | детерминированность манифеста прогона |

### CUDA / RTX 4090

На машине с CUDA Toolkit 12.x и RTX 4090 собирайте независимые CPU и GPU
каталоги: `CMAKE_CUDA_ARCHITECTURES=89` уже задан в preset. CUDA fast math
отключён намеренно, пока не доказан parity. Не используйте GPU-бинарник как
CPU-reference.

```sh
cmake --preset cpu-release
cmake --build --preset cpu-release
cmake --preset cuda-release
cmake --build --preset cuda-release
ctest --preset cuda-release -E '^mpi\.decomposition_parity$'
python3 tests/check_cpu_gpu_parity.py \
  --cpu build/cpu-release/mhd2d --gpu build/cuda-release/mhd2d \
  --config inputs/uniform_const.json --output-dir benchmarks/raw/cuda/parity-uniform
python3 tests/check_cpu_gpu_parity.py \
  --cpu build/cpu-release/mhd2d --gpu build/cuda-release/mhd2d \
  --config inputs/orszag_tang_uniform.json --output-dir benchmarks/raw/cuda/parity-orszag
```

Parity сопоставляет `rho/p`, нормы `div B` и счётчики fallback/floor. До его
прохождения измерять скорость нельзя. Multi-rank MPI/GPU, Nsight и performance
conclusions этим single-GPU gate не покрываются.

---

## 3. Сценарии нового решателя

```sh
./build/release/mhd2d inputs/orszag_tang.json      # вихрь Орзага–Танга, AMR 2 уровня
./build/release/mhd2d inputs/rotor.json            # вращающийся цилиндр
./build/release/mhd2d inputs/brio_wu.json          # ударная труба Брио–Ву
./build/release/mhd2d inputs/alfven_wave.json      # CP-альфвеновская волна
./build/release/mhd2d inputs/magnetic_loop.json    # перенос петли поля
./build/release/mhd2d inputs/mhd_blast.json        # МГД-взрыв (декартов)
./build/release/mhd2d inputs/uniform_const.json    # постоянное состояние
mpirun -np 4 ./build/release/mhd2d inputs/orszag_tang.json
```

Каждый прогон печатает в конце строки `ranges:`, `divb:`, `conservation:`,
`regrid_jump:` и `Evolve finished:` — по ним и работают автотесты.

Быстрый автономный драйвер (без AMReX), удобен для одномерных задач:

```sh
./build/release/mhd2d_verify briowu1d 400 mc rk2 gs 0.1 out.csv
./build/release/mhd2d_verify ot 128        # -> out_ot.csv
./build/release/mhd2d_verify loop 128 2.0 0.1
./build/release/mhd2d_verify dw1d 400 none euler bs 0.2 dw.csv
./build/release/mhd2d_verify alfven 64
```

Аблация задаётся конфигом, а не пересборкой: `scheme.limiter`
(`none` = первый порядок), `time.integrator` (`euler` | `rk2`),
`scheme.emf_averaging` (`balsara_spicer` | `gardiner_stone`).

Автономный драйвер принимает восьмым аргументом набор переменных
реконструкции — `prim` (рабочий) или `cons`; последний нужен только для
проверки гипотезы о происхождении перелёта (§ 4 в `docs/RP1_MONOTONICITY.md`).

Полная матрица аблации «одна деталь за раз» против независимого эталона:

```sh
./build/release/briowu_reference 6400 0.1 benchmarks/raw/rp1_limiters/kt_ref_6400.csv 0.4
python3 scripts/limiter_study.py --verify ./build/release/mhd2d_verify \
    --reference-csv benchmarks/raw/rp1_limiters/kt_ref_6400.csv --nx 400 --cfl 0.1 \
    --output benchmarks/summary/rp1_limiter_study_n400.json
```

---

## 4. Сценарии исторического решателя

Всегда через раннер: он клонирует неизменяемый источник, накладывает оверлей,
собирает, гоняет CTest, запускает решатель и пишет манифест с quality gate —
даже если прогон провалился.

```sh
python3 scripts/run_legacy_corrected.py \
    --source /Users/ivansamanov/Documents/MHD2D \
    --case rotor --mesh-backend structured \
    --artifact-dir benchmarks/raw/legacy_corrected/rotor_128 \
    --compiler g++-15 --jobs 8
```

Карты: `brio_wu`, `cp_alfven`, `magnetic_loop_athena`,
`magnetic_loop_legacy_scaled`, `rotor`, `orszag_tang`.

Нерегулярная (Netgen) сетка — стресс-тест, на котором видны ошибки, невидимые
на структурной:

```sh
python3 -m venv .venv-netgen
.venv-netgen/bin/pip install numpy netgen-mesher
python3 scripts/run_legacy_corrected.py --source <MHD2D> --case rotor \
    --mesh-backend netgen --netgen-python .venv-netgen/bin/python --maxh 0.012 \
    --artifact-dir benchmarks/raw/legacy_corrected/rotor_netgen
```

Ключи конфига, специфичные для исправленного профиля:

| ключ | значения | смысл |
|---|---|---|
| `cfl` | `> 0` | отрицательное/отсутствующее сохраняет исторический умолчание карты |
| `ctEnergyMode` | `conservative` \| `preserve_internal` | что сохранять при RT0-реконструкции (D-009) |
| `pressureFloor` | `0` \| `> 0` | `0` — неположительное `p` останавливает расчёт; `> 0` — пол с подсчётом событий |
| `cylindrical`, `gpu` | — | **отказывают** с кодом 2: эти пути не исправлены |

---

## 5. Замеры и профиль

```sh
python3 scripts/benchmark.py --executable ./build/release/mhd2d \
    --config inputs/orszag_tang_uniform.json --label "OT 128^2" \
    --repeats 5 --warmup 1 --output benchmarks/summary/timing_ot128.json

python3 scripts/scaling.py --executable ./build/release/mhd2d \
    --config inputs/orszag_tang_uniform.json --mode omp --counts 1,2,4,8

cmake --build --preset profile -j 8
./build/profile/mhd2d inputs/amr_conservation.json   # печатает TinyProfiler
```

Замеры на одной незакреплённой рабочей станции **диагностические**; как
результаты масштабирования их приводить нельзя (см. `docs/T09_TIMING.md`).

---

## 6. Отчёт

```sh
sh scripts/regen_report_data.sh          # входные данные всех рисунков
python3 scripts/make_report_figures.py   # -> docs/figures/data/*.dat
python3 scripts/field_map.py benchmarks/raw/report_inputs/ot_128.csv rho \
    docs/figures/maps/ot_rho.png --upscale 4
cd docs && latexmk -pdf report.tex       # 19 страниц, 0 overfull
```

Генератор рисунков **отказывается** записывать в заголовок эфемерный путь
(`/tmp`): рисунок, происхождение которого нельзя воспроизвести, в отчёт не
попадёт молча.

---

## 7. Кластер и Ubuntu/RTX 4090

Перед первым запуском соберите preflight-отчёт на **compute node**, а не на
login node:

```sh
bash scripts/cluster/collect_system_info.sh --output /tmp/mhd-preflight.txt
```

В отчёте должны быть `nvidia-smi`, `nvcc`, компилятор, CMake, MPI, Slurm,
CPU/RAM/диск, загруженные modules и commit checkout. Пароли, токены и ключи
в него не попадают; пути/имя пользователя при необходимости можно редактировать
перед отправкой.

Доступа и параметров планировщика нет (решение D-005), поэтому скрипты в
`scripts/cluster/` — **шаблоны**, они не проверены на реальной очереди. Перед
первым запуском заполнить `--account`, `--partition` и модули окружения.

Для новой кампании используйте единый launcher, а не редактирование sbatch
файлов в дереве проекта:

```sh
cp scripts/cluster/campaign.env.example scripts/cluster/sites/k10.env
# заполните путь к durable storage, account/partition, modules и лимиты
scripts/cluster/submit_campaign.sh --site scripts/cluster/sites/k10.env \
    --legacy-source /path/to/MHD2D --dry-run
scripts/cluster/submit_campaign.sh --site scripts/cluster/sites/k10.env \
    --legacy-source /path/to/MHD2D --cuda-validation
```

Он ставит `amrex CPU validation -> strong/weak scaling` и независимый
`legacy_corrected` smoke run. Исходные деревья не модифицируются: сборки живут
в node-local temporary directory, а результаты — в новом каталоге campaign.
Сначала всегда проверьте `--dry-run`. Launcher принимает только чистый
checkout и job сверяет его commit перед сборкой, поэтому результаты не могут
тихо попасть от изменившегося после submit кода.

На удалённом Ubuntu-хосте (включая RTX 4090):

```sh
scripts/cluster/run_ubuntu4090.sh --legacy-source /path/to/MHD2D \
    --artifact-root /data/mhd-artifacts
```

`--cuda-validation` создаёт отдельные CPU-reference и CUDA builds, запускает
CUDA CTest (кроме multi-rank MPI parity) и сохраняет логи двух CPU/GPU
parity-case. Это single-GPU correctness gate, а не GPU benchmark.

Что нужно получить от владельца до запуска: планировщик и его версия, account
и partition, лимиты по времени и памяти, спецификация узла (ядра, сокеты,
память), наличие и модель GPU, доступные модули компилятора и MPI.
