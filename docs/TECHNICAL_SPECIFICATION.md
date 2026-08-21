# Детальное численное ТЗ на завершение и исследование MHD-решателей

Статус: проект ТЗ, подготовленный по состоянию репозитория на 2026-08-21 и
уточнённый после аудита legacy-проекта и ВКР владельца.

Системный план программы из двух продуктов приведён в
[`PROGRAM_PLAN.md`](PROGRAM_PLAN.md). Настоящий документ детализирует численную
верификацию нового `mhd2d_amrex` и его сопоставление с `mhd2d_legacy`.

## 1. Цель работы

Поставить два рабочих решателя и общий контур проверки:

- `mhd2d_legacy` — воспроизводимую, исправленную и оптимизированную версию кода
  ВКР на неструктурированных сетках;
- `mhd2d_amrex` — новый решатель с доказанной корректностью, консервативным AMR
  и CPU/GPU масштабированием;
- единые постановки, метрики и benchmark protocol.

Для нового 2D/2.5D решателя идеальной магнитной гидродинамики требуется:

1. объяснить и количественно исследовать «немонотонность» в задаче Brio–Wu;
2. провести прямое сравнение прежней и новой численных схем;
3. доказать корректность однородных и AMR-расчётов;
4. измерить время, strong/weak scaling, эффективность и ускорение CPU-кода;
5. реализовать и проверить NVIDIA CUDA backend через AMReX;
6. подготовить кластерные job scripts, сырые результаты, таблицы и итоговый
   научно-технический отчёт.

ТЗ ориентировано на автономную работу сильных кодовых агентов. Любой численный
вывод должен иметь воспроизводимый вход, сырой результат и скрипт агрегации.

## 2. Что уже есть в новом проекте

### 2.1. Архитектура

Основной код занимает около 3600 строк C++/Python/LaTeX и состоит из:

| Компонент | Назначение |
|---|---|
| `src/kernels/MhdState.H` | консервативные/примитивные переменные, УРС, fast speed |
| `src/kernels/Hlld.H` | HLLD Miyoshi–Kusano, поток и граневая ЭДС |
| `src/kernels/Reconstruction.H` | MUSCL, MC/minmod/van Leer |
| `src/kernels/CtUpdate.H` | Balsara–Spicer и упрощённая Gardiner–Stone corner EMF |
| `src/Problems.H` | Brio–Wu, Alfvén wave, Orszag–Tang, rotor, magnetic loop |
| `src/MhdAmr.H/.cpp` | `AmrCore`, FillPatch, AMR, RK2, plotfile/HDF5 |
| `tests/standalone_verify.cpp` | uniform-grid драйвер тех же kernel headers |
| `scripts/plot_verification.py` | существующая визуализация CSV |
| `docs/REPORT.md`, `docs/report.tex` | текущий отчёт и его PDF |

Состояние хранится следующим образом:

- `rho`, три компоненты импульса, полная энергия и `Bz` — cell-centered;
- `Bx`, `By` — face-centered;
- `Ez` — nodal;
- клеточные `Bx`, `By` синхронизируются полусуммой противоположных граней.

На гранях выполняется MUSCL-реконструкция примитивных переменных; нормальная
компонента `B` берётся из face-centered массива. HLLD даёт газовый поток и ЭДС.
CT обновляет `Bx/By` дискретным законом Фарадея. На всех AMR-уровнях используется
общий `dt`, то есть временного subcycling нет.

### 2.2. Текущая проверенная исходная точка

- Release-сборка GCC 15 + AMReX 25.01, MPI + OpenMP, double precision проходит.
- Локальный `mhd2d_verify briowu`: 488 шагов при CFL 0.4, около 0.26 s на
  исследованной машине, `max|divB| = 0`.
- Локальный основной `mhd2d` для Brio–Wu: 391 шаг при CFL 0.5, около 0.70 s
  end-to-end с инициализацией MPI и plotfile I/O. Это диагностические числа,
  **не benchmark**: выполнен один запуск без pinning и статистики.
- В текущем CSV Brio–Wu для MC наблюдаются `rho_min ≈ 0.1166` при правом
  начальном состоянии 0.125 и `p_min ≈ 0.0870` при правом состоянии 0.1.
  Эти значения указывают на undershoot, но сами по себе не отделяют физическую
  структуру решения от паразитной осцилляции.
- Git-ветка `main` имеет baseline-коммиты `b0d1278` и `4dc9ad1`; перед каждым
  work package требуется фиксировать актуальный commit и dirty-state.

## 3. Результаты аудита и исходные риски

### 3.1. Критические (P0)

1. **Нет conservative reflux для газовых величин на AMR.**
   `SyncEmfAcrossLevels()` согласует nodal `Ez`, но `flux_` разных уровней не
   проходит через `FluxRegister` и не заменяется area-averaged fine flux.
   `AverageDownAll()` не компенсирует flux mismatch. До исправления нельзя
   заявлять глобальную консервативность AMR.

2. **GPU отсутствует как исполняемый путь.**
   В `MhdAmr.cpp` вычисления, инициализация, ГУ, tagging, редукции и derived
   fields используют `LoopOnCpu`. CMake не настраивает приложение через
   `setup_target_for_cuda_compilation`. `AMREX_GPU_HOST_DEVICE` есть только в
   независимых заголовках и не запускает kernel.

3. **Нет воспроизводимого regression/benchmark harness.**
   CTest отсутствует; параметры uniform driver захардкожены; нет machine-readable
   метрик, provenance, strong/weak scaling scripts и статистики повторов.

4. **Legacy-код и внешний архив найдены, но exact run пока невоспроизводим.**
   Зафиксирован `MHD2D` commit
   `9d0f60ea8576fac5d6f28c4dec142236d76131d6`. Inputs, Netgen mesher и 60 GB
   результатов найдены на `/Volumes/Elements`, но лежат в иных dirty snapshots
   без run manifest. Локальная сборка baseline на AppleClang останавливается
   на `omp.h`. См. [аудит](LEGACY_BASELINE.md) и
   [external manifest](EXTERNAL_DATA_MANIFEST.md).

### 3.2. Высокие (P1)

1. HLLD не имеет positivity-preserving fallback на HLLE/Rusanov при
   вырожденных/нефизичных intermediate states. Floors в `cons_to_prim` скрывают
   проблему, но не исправляют консервативное состояние.
2. В горячем AMReX-пути используются локальные `double`, хотя интерфейс проекта
   декларирует `mhd::Real`. Single precision build не является поддержанным.
3. Полная энергия и cell-centered магнитная энергия после CT могут быть
   несогласованы. Политика correction отсутствует и должна быть исследована.
4. Формула corner EMF названа Gardiner–Stone, но является сокращённой формой.
   Требуется проверить одномерную согласованность, устойчивость и отличие от
   полноценного upwind CT.
5. Нет тестов сохранения интегральной массы/импульса/энергии и AMR/uniform
   эквивалентности.
6. В конце расчёта plotfile иногда записывается дважды: `plot_now` на последнем
   шаге и безусловный `WritePlotFile` после цикла. AMReX переименовывает первый
   каталог в `.old.*`, что и создаёт лишний I/O.

### 3.3. Средние (P2)

- Все примеры запрашивают HDF5, но default build компилируется без HDF5 и молча
  откатывается к native plotfile после warning.
- Нет checkpoint/restart, что рискованно для длинных кластерных запусков.
- Родительский `output_dir` создаёт только IO rank без явного barrier.
- В документации есть неподтверждённые формулировки о «полной
  верифицированности», GPU-готовности и разрешении фронтов 2–3 ячейки.
- Существующие AMReX и standalone постановки отличаются CFL и интерфейсом
  параметров, поэтому их результаты нельзя смешивать в одной таблице.
- Артефакты plotfile/build/results находятся в корне и пока не отделены от
  исходников правилами Git.

## 4. Уточнение требований заказчика

### 4.1. «Вместо SSP взять RK2; TVD Рунге–Кутта (Шу)»

Текущий код уже реализует

```text
U(1)   = U(n) + dt L(U(n))
U(n+1) = 1/2 U(n) + 1/2 [U(1) + dt L(U(1))]
```

Это explicit two-stage SSP-RK2/Heun, то есть TVD RK2 Шу–Ошера при условии,
что forward Euler spatial operator TVD при допустимом CFL. SSP — современное
название более общего свойства ранее называвшихся TVD time discretizations.

Требуемое действие: не заменять формулу, а:

1. назвать enum и документацию однозначно `ssprk2_shu_osher`;
2. оставить alias `rk2` для обратной совместимости;
3. добавить unit test коэффициентов Бутчера и regression order test;
4. сравнить Euler и SSP-RK2 в ablation;
5. не добавлять midpoint RK2 как default, поскольку он не имеет той же SSP
   гарантии для нелинейного TVD-оператора.

### 4.2. Что означает «немонотонность Brio–Wu»

Задача не обязана иметь монотонный профиль каждой переменной: её решение
содержит несколько семейств волн, составную волну и физические экстремумы.
Кроме того, идеальная МГД не строго гиперболична; для классических данных
Brio–Wu существуют regular и семейство non-regular решений. Поэтому нужно
исследовать три разных вопроса:

1. какую ветвь решения выбирает HLLD + данная дискретизация;
2. какие экстремумы присутствуют в согласованном reference solution;
3. какой excess variation/overshoot порождает сетка поверх эталона.

### 4.3. «Посчитать на кластере»

ТЗ включает scripts и протокол, но фактический расчёт возможен только после
предоставления доступа, scheduler/account/partition и характеристик узлов.
Таблицы нельзя заранее заполнять оценочными числами.

### 4.4. Работа Бисикало–Жилкина

По формулировке не удалось однозначно определить публикацию. Для темы
монотонности наиболее релевантна работа Дудорова–Жилкина–Кузнецова 1999 года
о явной квазимонотонной conservative TVD-схеме. Работы Бисикало–Жилкина в
основном применяют MHD-модели к аккреционным течениям. Заказчик должен дать
точное название или PDF; до этого обе линии литературы считать кандидатами.

## 5. Обязательные результаты проекта

В конце программы должны быть поставлены:

1. самостоятельный `mhd2d_legacy` с профилями `vkr`, `corrected`, `opt_cpu`,
   `opt_cuda`, mesh pipeline и parity tests;
2. самостоятельный `mhd2d_amrex` с uniform/AMR и CPU/MPI/OpenMP/CUDA;
3. корректный Git baseline и `.gitignore` для каждого source tree;
4. CTest suite для kernels, uniform grid, AMR, MPI и CUDA;
5. конфигурируемый verification/benchmark runner;
6. эталон и машинные метрики Brio–Wu;
7. harness сравнения legacy/new;
8. conservative AMR synchronization для gas + CT;
9. CPU timers и raw benchmark format;
10. strong/weak scaling scripts и таблицы;
11. CUDA builds, GPU-порты и CPU/GPU parity tests;
12. Slurm scripts для CPU, hybrid и GPU запусков;
13. итоговые CSV/JSON, графики, таблицы Markdown/LaTeX;
14. обновлённый научно-технический отчёт без неподтверждённых утверждений.

## 6. Детальные пакеты работ нового solver и сравнения

Legacy-пакеты L0–L3 и общая последовательность программы определены в
`PROGRAM_PLAN.md`. Ниже приведена детализация N0–N5/C0–C2 для текущего
репозитория.

### WP0. Воспроизводимая исходная точка

Задачи:

- использовать существующий baseline commit и создавать отдельный проверяемый
  commit для каждого связного work package;
- добавить `.gitignore` для `build*/`, `results/`, `plt_*`, `*.bp`, временных
  профилей и локальных бинарников;
- не удалять существующие пользовательские артефакты; при необходимости
  перечислить их в manifest;
- добавить `CMakePresets.json`: `cpu-debug`, `cpu-release`, `mpi-release`,
  `cuda-release`, `profile`;
- включить `enable_testing()` и зарегистрировать быстрые CTest;
- печатать в каждый run manifest: commit, dirty flag, compiler, flags, AMReX,
  MPI/OpenMP/CUDA versions, input SHA-256, hostname и UTC timestamp;
- валидировать JSON: `gamma>1`, положительные размеры/CFL, допустимые enum,
  парность periodic BC, совместимость output format со сборкой;
- добавить `benchmark_mode`, полностью отключающий plot/diag I/O, но не
  критические проверки finite/positivity;
- устранить двойную финальную запись plotfile.

Критерий приёмки:

- чистая configure/build/test последовательность воспроизводится из нового
  checkout;
- CTest падает при испорченном HLLD/RK/CT kernel;
- два запуска с одним input дают manifest с одинаковой численной конфигурацией.

### WP1. Unit- и regression-тесты численного ядра

Добавить tests для:

1. `prim_to_cons` ↔ `cons_to_prim` на случайных физических состояниях;
2. fast magnetosonic speed в гидро-, Alfvén- и вырожденных пределах;
3. HLLD consistency `F(q,q)=F_phys(q)`;
4. supersonic left/right branches HLLD;
5. rotational/contact discontinuity и `Bn→0`;
6. конечности всех flux components на параметрическом наборе states;
7. limiter symmetry, constant preservation и отсутствие новых scalar extrema;
8. corner EMF constant-state и одномерный предел;
9. CT: неизменность дискретного `div B` после случайного nodal `Ez` update;
10. SSP-RK2 second-order convergence на автономной гладкой ODE;
11. manufactured/advection smooth test для совместного порядка пространства и
    времени;
12. Alfvén convergence минимум на N = 16, 32, 64, 128.

Для негативных `rho/p`, NaN или HLLD degeneration:

- добавить явный validity predicate;
- считать число реконструкционных fallback и flux fallback;
- fallback минимум на HLLE с консервативным потоком;
- результаты fallback counter писать в diagnostics/manifest;
- acceptance: на стандартных гладких тестах fallback = 0, на стресс-тестах
  решение остаётся finite и положительным.

### WP2. Исследование Brio–Wu

#### WP2.1. Эталон

Зафиксировать один основной reference branch и один sensitivity branch:

- основной: классическое compound/non-regular решение, обычно выбираемое
  shock-capturing schemes;
- sensitivity: regular solution либо слегка non-coplanar perturbed problem,
  устраняющий часть вырождения.

Reference получить одним из способов в порядке предпочтения:

1. опубликованный exact Riemann solver, умеющий regular/non-regular waves;
2. независимый доверенный код на N ≥ 16384 с документированным solver;
3. grid-converged внутренний результат, только после сравнения с независимым
   источником.

Нельзя использовать исследуемый N=512 result как единственный эталон для себя.

#### WP2.2. Матрица расчётов

Uniform grid, `t=0.1`, `gamma=2`, `Bx=0.75`, одни и те же НУ/ГУ:

| Фактор | Значения |
|---|---|
| N по x | 128, 256, 512, 1024, 2048, при необходимости 4096 |
| y | минимальное число ячеек, подтверждающее 1D-инвариантность; отдельно 4 |
| CFL | 0.2, 0.4, 0.5 |
| limiter | none, minmod, MC, van Leer |
| time | Euler, SSP-RK2 Shu–Osher |
| corner EMF | Balsara–Spicer, current GS, улучшенный upwind CT при наличии |
| precision | double; single только как sensitivity после поддержки |

Главный набор — заранее выбранные N=256/512/1024, CFL=0.4, без перебора после
просмотра результата. Остальная матрица — sensitivity.

#### WP2.3. Метрики

Для `rho, u, p, By` вычислять:

- `L1`, `L2`, `Linf` против reference, спроецированного как cell averages;
- `TV(q_h)` и excess `max(0, TV(q_h)-TV(q_ref))`;
- локальные overshoot/undershoot относительно reference envelopes;
- число sign changes первой разности после исключения эталонных экстремумов;
- положение каждой волны и ошибка положения в `dx`;
- ширина скачка `N_10-90` в ячейках и физических единицах;
- peak/plateau error у составной волны и контакта;
- `min(rho)`, `min(p)`, NaN count, fallback counters;
- `max|divB|` и нормированную `dx max|divB| / max|B|`;
- wall time и cell-updates/s для compute-only запуска.

Фронты размечаются по reference один раз; разметка хранится в YAML/JSON.
Визуально утверждать «2–3 ячейки» можно только после `N_10-90` table.

#### WP2.4. Диагноз и возможные исправления

Проверять по одному фактору:

1. component-wise primitive limiting против characteristic limiting;
2. MC против minmod и shock flattening;
3. HLLD wave speed estimate и positivity fallback;
4. SSP CFL coefficient;
5. 1D-consistency corner EMF;
6. energy correction после CT;
7. влияние `ny`, BC и двухмерного пути на строго 1D solution.

Исправление принимается, если excess variation/overshoot уменьшается без
статистически значимого ухудшения фронтов и гладкой second-order convergence.

### WP3. Прямое сравнение старой и новой схем

#### WP3.1. Зафиксированная legacy-база

Исходная точка сравнения:

- каталог `/Users/ivansamanov/Documents/MHD2D`;
- remote `https://github.com/Ship-Vano/MHD2D`;
- clean `master` на commit
  `9d0f60ea8576fac5d6f28c4dec142236d76131d6`;
- ВКР `VKRB_Shamanov_FINITE.pdf`, SHA-256
  `bc5b2661fd65c8ac94f2a02e0236978515af8263da86c857d1a6efb9bea67645`;
- полный технический аудит: `docs/LEGACY_BASELINE.md`.
- внешний data/source manifest: `docs/EXTERNAL_DATA_MANIFEST.md`.

Перед запуском требуется сопоставить найденные external input/mesh/result с
конкретным commit и case. Оригинальные snapshots не модифицировать; build
fixes, adapters и neutral cases хранить отдельно в `mhd-amrex`.

#### WP3.2. Нормализация

Создать neutral case description и adapters для обоих кодов. Сверить:

- нормировку `B` и наличие `4π`;
- total energy definition;
- cell centers/face averages;
- координату исходного разрыва;
- BC и ghost policy;
- CFL definition (directional minimum или sum of spectral radii);
- time of sampling;
- precision и compiler optimizations.

Дополнительно учесть установленные различия:

- legacy — piecewise constant + forward Euler на треугольниках;
- new — MUSCL + SSP-RK2 на Cartesian grid;
- legacy Brio–Wu в текущем 2D-коде принудительно использует CFL `0.9`, тогда
  как ВКР использовала `0.1`, а текущий new JSON — `0.5`;
- legacy использует free-flow 2D boundaries, ВКР описывает frozen boundaries,
  new — fixed x и periodic y;
- legacy `computeDivergence()` не берёт модуль, поэтому его напечатанный
  максимум нельзя напрямую сравнивать с `max|divB|` нового кода.
- CPU legacy HLLD усредняет нормальный `B`, а CUDA-копия отличается от CPU
  оценкой правой скорости волны; сравнение требует CPU/GPU flux parity test.

#### WP3.3. Обязательные сравнения

1. **Equal spacing:** одинаковые `dx, dy` и `t_end`.
2. **Equal degrees of freedom:** отдельно, если старый grid topology иной.
3. **Equal cost:** подобрать N по median compute time.
4. **Equal error:** интерполировать кривую time-vs-error.

Набор задач:

- Brio–Wu — разрывы и осцилляции;
- Alfvén wave — порядок, диссипация, фазовая ошибка;
- Orszag–Tang — сложное 2D течение и reference slices/maps;
- rotor — сильное magnetic braking и устойчивость;
- magnetic loop — advection/diffusion магнитного поля;
- AMR conservation test — только после WP4.

Для новой схемы обязательно выполнить ablation:

| ID | Space | Time | EMF |
|---|---|---|---|
| N0 | piecewise constant | Euler | Balsara–Spicer |
| N1 | MUSCL selected limiter | Euler | Balsara–Spicer |
| N2 | MUSCL | SSP-RK2 | Balsara–Spicer |
| N3 | MUSCL | SSP-RK2 | selected upwind CT |

Такое разложение не позволяет приписать всё улучшение только AMReX или HLLD.

Сначала выполнить контролируемое 1D-сравнение legacy HLLD с N0 на едином
sampling/projection protocol и точно одинаковом output time. Затем сравнивать
полные 2D-коды, проецируя треугольные данные
area-weighted на общие x-bins. Историческую таблицу ошибок ВКР и её CPU/GPU
времена помечать `historical_unreproduced`, пока не получены raw logs.

#### WP3.4. Выходная таблица качества

```text
case, solver, commit, mode, N, dof, steps, dx, CFL,
L1, L2, Linf, front_width_cells, front_position_error_dx,
overshoot, undershoot, TV_excess, max_divB, min_rho, min_p,
compute_s, end_to_end_s, fallback_count
```

### WP4. Корректный AMR

Задачи:

1. добавить `FluxRegister` для `rho`, momentum, energy, `Bz` либо до update
   заменить coarse interface flux точным area-average fine flux;
2. не reflux cell-centered copies `Bx/By`; они определяются face CT;
3. сохранить EMF synchronization, но проверить знаки/геометрию на всех
   orientation coarse–fine interfaces;
4. добавить интегральные diagnostics с masking covered coarse cells;
5. проверить regrid interpolation и average-down без изменения conservation;
6. рассмотреть checkpoint/restart;
7. документировать, почему subcycling выключен; subcycling не добавлять до
   готовности time-integrated flux и EMF registers.

AMR regression cases:

- advection гладкой плотности через неподвижную refinement boundary;
- MHD wave через coarse–fine interface в обе стороны;
- static divergence-free field при нескольких regrid;
- fixed-grid two-level case для точного сравнения uniform fine grid;
- Orszag–Tang с AMR и uniform finest-equivalent reference.

Acceptance:

- conservation defect соответствует roundoff + физическому boundary flux;
- `div B` остаётся roundoff-scaled, без скачка при regrid;
- ошибка AMR не хуже uniform coarse и уменьшается к uniform fine;
- MPI decomposition не меняет решение сверх заданного floating-point tolerance.

### WP5. Инструментирование времени

Включить `AMReX_TINY_PROFILE=ON` в profile preset и добавить `BL_PROFILE` для:

- initialization;
- FillPatch cells/faces;
- physical BC;
- primitive/reconstruction + HLLD fluxes;
- corner EMF;
- flux/EMF synchronization;
- cell/CT update;
- average-down/reflux;
- ComputeDt/diagnostics;
- regrid;
- plot/checkpoint I/O.

Дополнительно писать один JSON/CSV record на measured run. Времена брать после
MPI barrier в начале и перед чтением wall time в конце; reported wall — max по
рангам. Для GPU учитывать асинхронность: timer вокруг всего `MFIter` region либо
явная stream synchronization только в benchmark instrumentation.

Два режима:

- `compute-only`: no plot/checkpoint, diagnostics вне timed region;
- `end-to-end`: production I/O и diagnostics включены.

Минимум 1 warm-up + 5 repeats; для коротких задач увеличить размер/число шагов,
чтобы run был не короче примерно 10 s на исследуемом ресурсе. Хранить все
повторы. Основной агрегат — median; разброс — min/max и MAD либо IQR.

### WP6. CPU parallel scaling

#### WP6.1. Strong scaling

Фиксированные задачи:

- uniform Orszag–Tang 2048² или больше;
- fixed two-level AMR с неизменным layout/coverage;
- production dynamic AMR отдельно.

Ресурсы, если позволяют узлы:

- MPI: 1, 2, 4, 8, 16, 32, 64, ... ranks;
- OpenMP на одном узле: 1, 2, 4, 8, ... threads;
- hybrid: 1/2/4 ranks per socket × threads до полного узла;
- повторить минимум на 1, 2, 4, 8 узлах.

Формулы:

```text
S(p) = median(T_ref) / median(T_p)
E(p) = S(p) / (p / p_ref)
```

`p` явно определяется как CPU cores для CPU-only experiment. SMT on/off —
отдельные серии. Фиксировать `OMP_PROC_BIND`, `OMP_PLACES` и launcher binding.

#### WP6.2. Weak scaling

Поддерживать примерно постоянные cells/core и одинаковое число steps. Domain
увеличивать так, чтобы aspect ratio менялся минимально. Для AMR фиксировать
coverage или использовать deterministic tagging snapshot.

```text
E_weak(p) = median(T_ref) / median(T_p)
```

#### WP6.3. Таблица

| nodes | ranks | threads/rank | cores | cells | levels | median compute s | IQR s | speedup | efficiency | cell-updates/s | comm % | regrid % | I/O s |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|

Дополнительно построить:

- time и speedup vs cores;
- parallel efficiency vs cores;
- time breakdown stacked plot;
- load imbalance `max/mean rank compute time`;
- AMR cells per rank distribution.

### WP7. CUDA/GPU порт

#### WP7.1. Сборка

- добавить option/preset `AMReX_GPU_BACKEND=CUDA`;
- перед добавлением AMReX включить CUDA language при требовании версии;
- после задания sources применить `setup_target_for_cuda_compilation(mhd2d)`;
- параметризовать `AMReX_CUDA_ARCH`;
- оставить один source path для CPU и GPU;
- в CI хотя бы compile-only CUDA job, если runner без GPU.

#### WP7.2. Порт циклов

Перевести на `ParallelFor`/`ParReduce`:

- InitLevelData;
- ErrorEst;
- face physical BC;
- SyncCellB;
- flux/reconstruction/HLLD;
- corner EMF;
- CT/gas update;
- ComputeDt и MaxDivB;
- derived plot fields.

`Problem`/configuration data, нужные device, представить trivially copyable POD
и device-callable functions/switch, без `std::function`. CPU boundary callback
заменить на GPU-capable boundary functor либо отделённый staging path с явно
измеренными transfers.

Использовать `mhd::Real` во всех device arrays. Проверить HLLD local arrays,
register pressure, nested lambdas, `std::fabs/sqrt` и device compilation.

#### WP7.3. GPU correctness

- CPU double и GPU double: одинаковые inputs и steps;
- нормы разности всех полей на каждом checkpoint;
- одинаковые conservation/divB trends;
- одинаковые fallback counters;
- compute-sanitizer на малом тесте;
- 1 GPU и multi-GPU MPI.

Tolerance задаётся до эксперимента и масштабируется с числом операций; bitwise
identity не требуется из-за порядка редукций.

#### WP7.4. GPU performance

- 1 MPI rank per GPU как основной layout;
- отдельный тест GPU-aware MPI on/off;
- strong scaling 1, 2, 4, 8, ... GPU;
- CPU-vs-GPU сравнивать как time-to-solution на одном полном CPU node против
  одного GPU и полного GPU node;
- исключить problem sizes, не насыщающие GPU, из главной таблицы, но сохранить
  их как latency results;
- Nsight Systems: kernel/communication/transfer timeline;
- Nsight Compute: occupancy, register pressure, memory throughput, roofline для
  HLLD/reconstruction и update kernels.

GPU speedup:

```text
S_gpu = median(T_best_full_CPU_node) / median(T_GPU_resource)
```

Не сравнивать один CPU core с целым GPU как основной маркетинговый результат.

### WP8. Кластерные расчёты

Создать шаблоны:

- `cluster/slurm/build_cpu.sh`;
- `cluster/slurm/build_cuda.sh`;
- `cluster/slurm/strong_cpu.sbatch`;
- `cluster/slurm/weak_cpu.sbatch`;
- `cluster/slurm/strong_gpu.sbatch`;
- `cluster/slurm/profile_gpu.sbatch`;
- `benchmarks/collect.py` и `benchmarks/summarize.py`.

Каждый job пишет отдельный run directory и manifest. Job arrays допустимы для
resource matrix, но повторы разных размеров должны иметь независимые IDs.

До запуска запросить:

- адрес/login method;
- Slurm/PBS/другой scheduler;
- account, partition/qos, wall limits;
- CPU model, cores/socket, NUMA, RAM;
- GPU model/count/node, compute capability;
- MPI, compiler, CUDA, HDF5 modules;
- GPU-aware MPI status;
- scratch/project paths и filesystem policy;
- разрешённые profiling tools.

### WP9. Итоговый отчёт

Отчёт должен разделять:

1. математическую постановку и нормировку;
2. численную схему и точную формулу SSP-RK2;
3. доказанные дискретные инварианты;
4. Brio–Wu: физическая структура, non-uniqueness, паразитные осцилляции;
5. legacy/new quality comparison;
6. convergence tables;
7. AMR conservation;
8. CPU strong/weak scaling;
9. GPU port, parity и speedup;
10. ограничения и отрицательные результаты.

Все рисунки генерируются одной командой из raw data. Таблицы в Markdown и
LaTeX должны строиться из одного summary CSV, чтобы числа не расходились.

## 7. Структура артефактов

Целевая структура:

```text
tests/
  unit/
  regression/
  reference/
benchmarks/
  cases/
  raw/                 # можно хранить вне Git, manifest обязателен
  summary/
  figures/
  scripts/
cluster/slurm/
docs/
  PROGRAM_PLAN.md
  TECHNICAL_SPECIFICATION.md
  REPORT.md
```

Идентификатор run:

```text
<utc>-<case>-<solver>-<commit>-n<nodes>-r<ranks>-t<threads>-g<gpus>-rep<k>
```

Raw record минимум:

```json
{
  "case": "orszag_tang_uniform",
  "solver": "new",
  "commit": "...",
  "input_sha256": "...",
  "nodes": 1,
  "ranks": 8,
  "threads_per_rank": 4,
  "gpus": 0,
  "steps": 100,
  "cells_updated": 0,
  "timing": {
    "compute_s": 0.0,
    "communication_s": 0.0,
    "regrid_s": 0.0,
    "io_s": 0.0,
    "end_to_end_s": 0.0
  },
  "quality": {
    "max_divB": 0.0,
    "min_rho": 0.0,
    "min_p": 0.0,
    "mass_defect": 0.0,
    "energy_defect": 0.0,
    "fallback_count": 0
  }
}
```

## 8. Gates нового solver и порядок выполнения

Работы выполнять в таком порядке:

| Gate | Содержание | Зависимость |
|---|---|---|
| G0 | Git baseline, presets, CTest, manifest, benchmark mode | — |
| G1 | kernel correctness + SSP-RK2 tests + positivity diagnostics | G0 |
| G2 | Brio–Wu reference/metrics и диагноз немонотонности | G1 |
| G3 | gas reflux + AMR conservation regression | G1 |
| G4 | legacy adapter и quality comparison | legacy inputs + подтверждённый commit, G2 |
| G5 | CPU profiling + local scaling harness | G0, желательно G3 |
| G6 | CUDA port + CPU/GPU parity | G1, G3 |
| G7 | cluster CPU/GPU campaigns | доступ к кластеру, G5, G6 |
| G8 | итоговые таблицы и отчёт | G2–G7 |

Запрещено переходить к публичным performance claims до G1 и к AMR performance
claims до G3. GPU speedup не публиковать до parity tests G6.

## 9. Definition of Done нового solver и общего эксперимента

Проект считается завершённым, когда одновременно выполнено:

- все обязательные CTest проходят CPU Release, MPI и CUDA;
- Brio–Wu metrics воспроизводятся из reference data, а причина каждого
  существенного экстремума классифицирована;
- legacy/new comparison содержит equal-spacing, equal-cost и equal-error views;
- гладкий тест показывает ожидаемый порядок, разрывной — convergence trend;
- AMR conservation defect и `divB` удовлетворяют заранее заданным tolerances;
- benchmark содержит raw repeats и машинный manifest;
- strong/weak scaling таблицы заполнены реальными кластерными данными;
- GPU port не содержит `LoopOnCpu` в timed hot path и проходит parity tests;
- CPU/GPU speedup определён относительно честного resource baseline;
- REPORT генерируется из тех же summary tables;
- известные ограничения перечислены явно, без слов «готово» для непроверенных
  возможностей.

## 10. Литературная база

Обязательное ядро:

1. M. Brio, C. C. Wu, *An upwind differencing scheme for the equations of
   ideal magnetohydrodynamics*, JCP 75 (1988),
   [DOI](https://doi.org/10.1016/0021-9991(88)90120-9).
2. T. Miyoshi, K. Kusano, *A multi-state HLL approximate Riemann solver for
   ideal MHD*, JCP 208 (2005),
   [DOI](https://doi.org/10.1016/j.jcp.2005.02.017).
3. E. N. Avdeeva, V. V. Lukin, *Divergence-free finite-difference method for
   2D ideal MHD*, JPCS 1336 (2019),
   [DOI](https://doi.org/10.1088/1742-6596/1336/1/012026).
4. C.-W. Shu, S. Osher, *Efficient implementation of essentially
   non-oscillatory shock-capturing schemes*, JCP 77 (1988),
   [DOI](https://doi.org/10.1016/0021-9991(88)90177-5).
5. S. Gottlieb, C.-W. Shu, *Total variation diminishing Runge–Kutta schemes*,
   Math. Comp. 67 (1998),
   [DOI](https://doi.org/10.1090/S0025-5718-98-00913-2).
6. S. Gottlieb, C.-W. Shu, E. Tadmor, *Strong Stability-Preserving High-Order
   Time Discretization Methods*, SIAM Review 43 (2001),
   [DOI](https://doi.org/10.1137/S003614450036757X).
7. G. Tóth, *The div B = 0 constraint in shock-capturing MHD codes*, JCP 161
   (2000), [DOI](https://doi.org/10.1006/jcph.2000.6519).
8. D. S. Balsara, D. S. Spicer, *A staggered mesh algorithm using high order
   Godunov fluxes to ensure solenoidal magnetic fields*, JCP 149 (1999),
   [DOI](https://doi.org/10.1006/jcph.1998.6153).
9. T. Gardiner, J. Stone, *An unsplit Godunov method for ideal MHD via
   constrained transport*, JCP 205 (2005),
   [DOI](https://doi.org/10.1016/j.jcp.2004.11.016).
10. D. S. Balsara, *Divergence-free adaptive mesh refinement for MHD*, JCP 174
    (2001), [DOI](https://doi.org/10.1006/jcph.2001.6917).
11. K. Takahashi, S. Yamada, *Regular and non-regular solutions of the Riemann
    problem in ideal MHD*, [arXiv:1210.5584](https://arxiv.org/abs/1210.5584).
12. А. Е. Дудоров, А. Г. Жилкин, О. А. Кузнецов, *Квазимонотонная разностная
    схема повышенного порядка точности для уравнений магнитной гидродинамики*,
    Матем. моделирование 11:1 (1999),
    [Math-Net](https://www.mathnet.ru/links/1dcb22cf86d706250aa1e3b59db9e3ca/mm1056.pdf).
13. И. П. Шаманов, *Разработка и исследование программного комплекса для
    решения задач идеальной магнитной гидродинамики*, ВКР бакалавра, МГТУ
    им. Н. Э. Баумана, 2025. Проверенный локальный PDF и checksum приведены в
    `docs/LEGACY_BASELINE.md`.

Инфраструктура и производительность:

14. [AMReX: GPU support](https://amrex-codes.github.io/amrex/docs_html/GPU.html).
15. [AMReX: AmrCore and FluxRegister](https://amrex-codes.github.io/amrex/docs_html/AmrCore.html).
16. [AMReX profiling tools](https://amrex-codes.github.io/amrex/docs_html/AMReX_Profiling_Tools.html).
17. [NVIDIA CUDA C++ Best Practices: timing and metrics](https://docs.nvidia.com/cuda/cuda-c-best-practices-guide/).
18. [NVIDIA Nsight Systems User Guide](https://docs.nvidia.com/nsight-systems/UserGuide/).
19. [NVIDIA Nsight Compute Profiling Guide](https://docs.nvidia.com/nsight-compute/ProfilingGuide/).

Отдельно запросить у заказчика точную работу Бисикало–Жилкина. ВКР и исходный
репозиторий старой схемы уже идентифицированы и внесены в baseline.

## 11. Данные, которые нужно получить от владельца сейчас

1. Подтверждение legacy commit `9d0f60e…` и точное соответствие найденных на
   `Elements` input/mesh/result конкретным запускам ВКР.
2. Точное название или файл статьи Бисикало–Жилкина.
3. Какой именно график/переменная названы заказчиком «немонотонными».
4. Требуемая ветвь/эталон Brio–Wu, если она задана руководителем.
5. Доступ и характеристики кластера.
6. Срок, формат отчёта и обязательные размеры таблиц/рисунков.
7. Нужно ли сохранять обратную совместимость существующих JSON и plotfile.

До получения пунктов 1–5 агенты могут полностью выполнить G0, G1, большую
часть G2/G3/G5 и подготовить G4/G6/G7, но не должны выдавать исторические
таблицы ВКР за воспроизведённые legacy-результаты или имитировать cluster data.
