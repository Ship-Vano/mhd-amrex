# mhd2d — 2D-решатель идеальной МГД на AMReX с constrained transport

Решатель двумерных уравнений идеальной магнитной гидродинамики по схеме
Авдеевой–Лукина (2019) на блочно-структурированной адаптивной сетке AMReX:

* газовые величины (ρ, ρ**v**, e) и Bz — в центрах ячеек;
* нормальные компоненты **B** (Bx, By) — на гранях (staggered);
* HLLD-решатель Римана (Miyoshi & Kusano, 2005) даёт потоки газовых величин
  **и** ЭДС на гранях;
* закон Фарадея аппроксимируется напрямую через теорему Стокса по узловым
  ЭДС; на проверенных однородных сетках дискретная `div B` остаётся малой;
* интегратор SSP-RK2 (Хойна), CFL ≈ 0.5; MUSCL-реконструкция (MC/minmod);
* AMR (`amrex::AmrCore`) c бездивергентной интерполяцией граней
  (`face_divfree_interp`) и синхронизацией ЭДС между уровнями; gas reflux и
  общая AMR-консервативность ещё должны быть реализованы и проверены;
* параметры запуска — JSON-файл (без перекомпиляции);
* вывод — нативные plotfile'ы AMReX или HDF5, оба открываются в ParaView.

Подробное описание схемы, тестов и результатов — в `docs/REPORT.md`.
План завершения, критерии приёмки и протоколы сравнения/benchmark — в
`docs/TECHNICAL_SPECIFICATION.md`; аудит прежнего решателя и ВКР — в
`docs/LEGACY_BASELINE.md`, внешний архив — в
`docs/EXTERNAL_DATA_MANIFEST.md`; правила работы AI-агентов — в `AGENTS.md`.

## Структура проекта

```
mhd-amrex/
├── CMakeLists.txt           # сборка (Ubuntu/macOS; MPI/OpenMP/HDF5 — опции)
├── src/
│   ├── kernels/             # вычислительные ядра, НЕ зависят от AMReX
│   │   ├── MhdState.H       #   переменные, УРС, скорости волн
│   │   ├── Hlld.H           #   HLLD-поток (+ЭДС)
│   │   ├── Reconstruction.H #   MUSCL-лимитеры
│   │   └── CtUpdate.H       #   узловые ЭДС, constrained transport
│   ├── MhdAmr.H/.cpp        # AmrCore: шаг по времени, AMR, FillPatch, вывод
│   ├── Problems.H           # каталог начальных условий (через вект. потенциал)
│   ├── Config.H/.cpp        # JSON-конфигурация
│   └── main.cpp
├── inputs/                  # примеры JSON-конфигураций
├── tests/standalone_verify.cpp  # автономный драйвер на тех же ядрах
├── scripts/plot_verification.py # графики и метрики
└── docs/REPORT.md           # отчёт: схема, тесты, эталоны, графики
```

## Сборка

### Ubuntu

```bash
sudo apt install build-essential cmake git libopenmpi-dev   # libhdf5-openmpi-dev для HDF5
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
      -DMHD_MPI=ON -DMHD_OPENMP=ON -DMHD_HDF5=OFF
cmake --build build -j $(nproc)
```

AMReX (тег `25.01`) и nlohmann/json скачиваются и собираются автоматически
(FetchContent). Для системного AMReX: `-DMHD_FETCH_AMREX=OFF` (нужна сборка
AMReX с `AMReX_SPACEDIM=2`).

### macOS

```bash
brew install cmake libomp open-mpi          # + hdf5-mpi для HDF5
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j 8
```

OpenMP у Apple clang подхватывается через `libomp` (CMake `find_package(OpenMP)`).
Если не находится: `-DOpenMP_ROOT=$(brew --prefix libomp)`.

### HDF5-вывод

```bash
cmake -S . -B build -DMHD_HDF5=ON
```

AMReX соберётся с `AMReX_HDF5=ON` и при `"format": "hdf5"` в JSON будет писать
`WriteMultiLevelPlotfileHDF5` (один `.h5` на снимок — быстрая коллективная
запись через MPI-IO). При `"format": "native"` пишутся обычные plotfile'ы.
**ParaView (≥ 5.7) открывает нативные plotfile'ы AMReX напрямую** (тип
*AMReX/BoxLib Grid Reader*) — поэтому HDF5 не обязателен для визуализации.

## Запуск

```bash
# последовательный
./build/mhd2d inputs/orszag_tang.json

# MPI (8 рангов) + OpenMP (4 потока на ранг)
OMP_NUM_THREADS=4 mpirun -np 8 ./build/mhd2d inputs/orszag_tang.json
```

Примеры конфигураций:

| Файл | Задача | Сетка | AMR | ГУ |
|---|---|---|---|---|
| `inputs/orszag_tang.json` | вихрь Орзага–Танга | 256² | 2 уровня | периодические |
| `inputs/brio_wu.json` | ударная труба Брио–Ву | 512×4 | нет | dirichlet («исторические») |
| `inputs/alfven_wave.json` | альфвеновская волна, α=30° | 64×110 | нет | периодические |
| `inputs/rotor.json` | вращающийся цилиндр | 200² | 1 уровень | outflow |

Все поля JSON описаны в `src/Config.H`; ключевые:

```jsonc
{
  "problem": "orszag_tang",                  // имя задачи из src/Problems.H
  "problem_params": { "B0": 0.282 },         // параметры конкретной задачи
  "geometry": { "prob_lo": [0,0], "prob_hi": [1,1], "n_cell": [256,256] },
  "amr":  { "max_level": 2, "refine_grad_rho": 0.2, "refine_current": 0.5 },
  "bc":   { "x_lo": "periodic|outflow|reflect|dirichlet", ... },
  "time": { "cfl": 0.5, "t_max": 0.5, "integrator": "rk2|euler" },
  "scheme": { "gamma": 1.667, "limiter": "mc|minmod|vanleer|none",
              "emf_averaging": "gardiner_stone|balsara_spicer" },
  "output": { "plot_dt": 0.1, "prefix": "plt_ot", "format": "hdf5|native" }
}
```

## Быстрая верификация без AMReX

Ядра схемы можно проверить автономным драйвером (собирается одним `g++`):

```bash
g++ -O3 -std=c++17 -Isrc/kernels tests/standalone_verify.cpp -o verify
./verify briowu          # Брио–Ву, 512 ячеек
./verify ot 192          # Орзаг–Танг 192²
./verify alfven 32       # альфвеновская волна, L1/L2-ошибки
./verify rotor 128
python3 scripts/plot_verification.py ot out_ot.csv   # графики + max|divB|
```

## Просмотр результатов в ParaView

1. *File → Open* → каталог `plt_ot00610/` (нативный) → reader **AMReX/BoxLib
   Grid Reader**; либо файл `plt_ot00610.h5` (HDF5-сборка).
2. Отметить нужные поля (`rho`, `p`, `Bx`, `By`, `divB`), *Apply*.
3. Для AMR-иерархии ParaView показывает уровни как vtkOverlappingAMR.

## Литература

* E.N. Avdeeva, V.V. Lukin. *Divergence-free finite-difference method for 2D
  ideal MHD*, J. Phys.: Conf. Ser. 1336 (2019) 012026.
* T. Miyoshi, K. Kusano. *A multi-state HLL approximate Riemann solver for
  ideal MHD*, JCP 208 (2005) 315.
* D. Balsara, D. Spicer. JCP 149 (1999) 270 — staggered CT с ЭДС из потоков.
* T. Gardiner, J. Stone. JCP 205 (2005) 509 — усреднение ЭДС, CTU+CT.
* G. Tóth. JCP 161 (2000) 605 — обзор div B-методов, постановки тестов.
* AMReX: https://amrex-codes.github.io/amrex/docs_html/ (AmrCore, FillPatch,
  `face_divfree_interp`, plotfile/HDF5).
* ВКРБ Шаманова И.П. (2025) — постановки тестов, «исторические» ГУ и
  невоспроизведённые legacy-результаты; provenance приведён в
  `docs/LEGACY_BASELINE.md`.
