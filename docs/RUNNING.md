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
ctest --preset release          # 21 тест, ~15 c
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
| `cuda-release` | объявлен, **на этой машине не собирается** (нет GPU NVIDIA) |

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
| `config.*` | строгая схема конфигурации, режимы аблации реально влияют |
| `manifest.repeatable` | детерминированность манифеста прогона |

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

## 7. Кластер

Доступа и параметров планировщика нет (решение D-005), поэтому скрипты в
`scripts/cluster/` — **шаблоны**, они не проверены на реальной очереди. Перед
первым запуском заполнить `--account`, `--partition` и модули окружения.

```sh
sbatch scripts/cluster/mhd2d_strong_scaling.sbatch
sbatch scripts/cluster/legacy_corrected_case.sbatch
```

Что нужно получить от владельца до запуска: планировщик и его версия, account
и partition, лимиты по времени и памяти, спецификация узла (ядра, сокеты,
память), наличие и модель GPU, доступные модули компилятора и MPI.
