# MHD-DUAL: план работ и системное техническое задание

Статус: базовый проект ТЗ, уточнён 2026-08-26.

Продукты: `mhd2d_legacy`, `mhd2d_amrex`.

Детальное численное ТЗ: [`TECHNICAL_SPECIFICATION.md`](TECHNICAL_SPECIFICATION.md).

Аудит старого кода: [`LEGACY_BASELINE.md`](LEGACY_BASELINE.md).

Фазовый план для Terra Ultra:
[`TERRA_EXECUTION_PLAN.md`](TERRA_EXECUTION_PLAN.md).

## 1. Цель и результат программы

Нужно поставить не один экспериментальный код, а два независимо собираемых и
проверяемых продукта с общим контуром верификации:

1. **`mhd2d_legacy`** — восстановленная, отлаженная и оптимизированная версия
   неструктурированного решателя из ВКР. Она служит историческим baseline и
   самостоятельной программой для тестов.
2. **`mhd2d_amrex`** — новый 2D/2.5D решатель идеальной МГД: HLLD + MUSCL + CT
   + SSP/TVD RK2 Шу–Ошера, AMR, MPI/OpenMP и CUDA.
3. **Общий V&V/benchmark-контур** — единые постановки, конвертеры, метрики
   качества, timing/scaling protocol и генерация итоговых таблиц.

Готовность определяется воспроизводимой сборкой, автоматическими тестами,
численными критериями и сырыми измерениями, а не только графиками решения.

### 1.1. Исследовательская гипотеза

Основная гипотеза H1: для выбранного класса 2D ideal-MHD задач реализация на
AMReX с MUSCL, SSP-RK2, CT и корректным AMR способна дать более выгодное
соотношение «ошибка — время — память», чем legacy finite-volume solver первого
порядка на неструктурированной сетке.

H1 не считается истинной заранее. Она проверяется в режимах equal resolution,
equal cost и equal error. Допустимым научным результатом является также
частичное или полное опровержение: например, AMR выигрывает на локализованных
структурах, но проигрывает на гладкой глобальной волне из-за overhead. В отчёте
запрещено заменять такую проверку фразой «AMReX быстрее и точнее» без данных.

### 1.2. Ближайший результат для научного руководителя

Первый пользовательский deliverable — исправленная сокращённая
расчётно-пояснительная записка НИРС для МГТУ им. Н. Э. Баумана. Она должна
честно описывать состояние «solver работает, доказательная база ещё строится»,
не заявлять гарантированную бездивергентность/AMR-консервативность/GPU-порт до
соответствующих gates и отделять измеренные результаты от плана продолжения.

## 2. Границы проекта

В scope входят восстановление ВКР и mesh pipeline, исправление известных
дефектов legacy, CPU/OpenMP/CUDA-оптимизация, завершение AMR и GPU в новом коде,
прямое сравнение разрешающей способности, замеры на workstation/кластере и
воспроизводимый отчёт.

Без отдельного change request не входят resistive/Hall MHD, общий 3D solver,
MPI-порт legacy и принципиальная замена HLLD. MPI и масштабируемый multi-GPU —
обязанность AMReX-продукта. Архивные VTU/plotfile размером в десятки гигабайт
не включаются в Git: сохраняются manifest, checksum и компактные срезы.

Полная реализация discontinuous Galerkin (DG) также не входит в текущую фазу.
Входит архитектурный задел: границы модулей, контракты состояния/потока/
граничных условий/интегратора, ADR и тесты, не позволяющие FV-реализации
монопольно определять всю кодовую базу.

Ключевые правила:

- numerical bug fix и performance optimization выполняются разными change sets;
- legacy остаётся независимой неструктурированной реализацией, а не AMReX-клоном;
- сравниваются одинаковые физические задачи, но не внутренние массивы программ;
- каждый результат связан с commit, dirty-state, config, dependencies, hardware
  и сырым логом;
- ускорение принимается только вместе с контролем изменения численного ответа.

## 3. Целевые продукты

### 3.1. `mhd2d_legacy`

Legacy-продукт обязан иметь четыре явно различимых профиля:

| Профиль | Назначение | Правило эквивалентности |
|---|---|---|
| `legacy_vkr` | Воспроизведение алгоритма ВКР | Математика заморожена; разрешены build/I/O/diagnostics |
| `legacy_corrected` | Верифицированная исправленная схема | Каждая numerical delta записана в bug-fix ledger |
| `legacy_opt_cpu` | Оптимизированный serial/OpenMP | Parity с выбранным CPU reference |
| `legacy_opt_cuda` | Оптимизированный CUDA | Parity с выбранным corrected reference |

`legacy_vkr` сохраняет исторический метод: finite volume на неструктурированной
треугольной сетке, первый порядок, forward Euler, HLLD, edge-normal `B`, nodal
EMF и Raviart–Thomas-реконструкцию магнитного поля. CFL и остальные параметры
задаёт версионируемый case manifest, а не скрытая константа.

Различие CPU/CUDA в правой скорости HLLD, signed-оценка `div B` и другие
обнаруженные дефекты исправляются в `legacy_corrected`; исторический профиль и
опубликованные результаты задним числом не изменяются.

### 3.2. `mhd2d_amrex`

Целевая конфигурация:

- Cartesian finite volume на AMReX;
- cell-centered газовые величины и `Bz`, face-centered `Bx/By`, nodal `Ez`;
- HLLD, MUSCL с конфигурируемым limiter, constrained transport;
- двухстадийный SSP-RK2/Heun — TVD RK2 Шу–Ошера;
- AMR с flux register/reflux и согласованием staggered magnetic field;
- MPI/OpenMP CPU backend и AMReX CUDA backend;
- диагностируемые positivity fallback, floors и invalid states.

Требование «вместо SSP взять RK2/TVD Шу» не означает замену текущей формулы:
реализованный SSP-RK2 уже является TVD RK2 Шу–Ошера. Нужны тест порядка и
фиксация алгоритма, а не midpoint RK2.

### 3.3. Общий validation/benchmark contract

Это не третий solver. Контур предоставляет canonical case schema, адаптеры
нативных конфигураций, mesh manifests, конвертацию VTU/AMReX plotfile,
вычисление ошибок/TV/front metrics/балансов, benchmark runner и путь данных
`raw -> summary -> figures -> report`. Его schema не зависит от типов AMReX.

### 3.4. Будущее расширение — DG для 2D MHD

DG планируется как новый spatial discretization того же физического ядра, а не
как набор `if (dg)` внутри FV/AMR-кода. Целевая декомпозиция:

```text
physics (EOS, state, flux, waves)
        ↓
spatial operator contract ── FV operator
        │                    DG operator (future)
        ↓
time integrator (Euler / SSP-RK)
        ↓
mesh/backend services (AMReX hierarchy, parallel execution)
        ↓
diagnostics / I/O / V&V
```

Повторно используются только математически общие части: state layout,
primitive/conservative conversion, physical flux, wave-speed estimates,
boundary-state policy, time-step orchestration, diagnostics и case schema.
FV reconstruction/CT update и будущие DG basis/quadrature/volume/surface
operators остаются раздельными. Для divergence control DG требуется отдельное
численное решение; текущий staggered CT нельзя объявлять автоматически
пригодным для DG.

## 4. Функциональные требования

### 4.1. Общие

| ID | Требование |
|---|---|
| SYS-001 | Обе программы собираются из clean checkout документированной командой. |
| SYS-002 | Запуск принимает явный config и создаёт run manifest. |
| SYS-003 | Manifest содержит commit/dirty-state, compiler/flags, зависимости, backend, hardware, ресурсы и физические параметры. |
| SYS-004 | Некорректный config даёт понятную ошибку и ненулевой exit code. |
| SYS-005 | Output содержит физическое время и однозначные имена/единицы полей. |
| SYS-006 | Diagnostic, plot и benchmark I/O управляются независимо. |
| SYS-007 | Case, mesh и reference имеют version/checksum. |
| SYS-008 | Научные утверждения в отчёте трассируются до requirement, test, raw artifact и commit. |
| SYS-009 | Physics kernels не зависят от конкретного spatial discretization и orchestration AMReX, кроме явно документированных adapters. |

### 4.2. Legacy

| ID | Требование |
|---|---|
| LEG-001 | Source snapshot ВКР сохраняется immutable tag. |
| LEG-002 | CMake и runtime не зависят от абсолютных путей пользователя. |
| LEG-003 | Netgen mesh pipeline воспроизводим; эталонные сетки имеют manifest/checksum. |
| LEG-004 | Brio–Wu, Alfvén, magnetic loop и Orszag–Tang выбираются config без перекомпиляции. |
| LEG-005 | Backends `serial/openmp/cuda` выбираются явно. |
| LEG-006 | Corrected CPU/CUDA используют согласованные HLLD-формулы. |
| LEG-007 | CUDA не выделяет и не копирует всю геометрию на каждом flux call. |
| LEG-008 | Диагностика выдаёт `max(abs(div B))`, scaled norm, min `rho/p`, NaN/Inf и интегральные балансы. |
| LEG-009 | OpenMP покрывает подтверждённые профилировщиком горячие циклы. |
| LEG-010 | `legacy_vkr` воспроизводит утверждённые таблицы ВКР либо даёт доказуемое объяснение расхождения. |
| LEG-011 | Каждый optimized backend проходит parity suite. |

### 4.3. Новый AMReX-код

| ID | Требование |
|---|---|
| NEW-001 | Uniform CPU solver имеет CTest для kernels и canonical cases. |
| NEW-002 | HLLD имеет устойчивый HLL/HLLE fallback для недопустимых промежуточных состояний. |
| NEW-003 | Floors/fallback имеют раздельные счётчики в diagnostics/manifest. |
| NEW-004 | Limiter, spatial order и ablation modes задаются config. |
| NEW-005 | SSP-RK2 общий для CPU/GPU и проверен temporal convergence test. |
| NEW-006 | Face flux и CT EMF согласованы в одномерном пределе. |
| NEW-007 | AMR реализует gas reflux и корректный sync/prolongation face `B`. |
| NEW-008 | Regrid проходит positivity, `div B` и conservation tests. |
| NEW-009 | Hot-path MultiFab loops используют `ParallelFor/ParReduce` и device-safe POD-параметры. |
| NEW-010 | CUDA build с явной architecture проходит тот же regression suite. |
| NEW-011 | Timing разделяет init, compute, communication/regrid, diagnostics, I/O и end-to-end. |
| NEW-012 | Один физический момент не выводится дважды без явного запроса. |

### 4.4. Прямое сопоставление

| ID | Требование |
|---|---|
| CMP-001 | Одинаковы НУ, ГУ, `gamma`, нормировка `B`, CFL policy, `t_end` и времена срезов. |
| CMP-002 | Выполняются режимы equal resolution, equal cost и equal error. |
| CMP-003 | Для разных сеток определена физическая мера `h`; числа ячеек/элементов публикуются отдельно. |
| CMP-004 | Для Brio–Wu заранее выбрана regular либо compound reference branch. |
| CMP-005 | Метрики: `L1/L2/Linf`, excess TV, overshoot/undershoot, паразитные экстремумы, положение/толщина фронтов. |
| CMP-006 | Ablation нового кода: first order; MUSCL; MUSCL+SSP-RK2; полный CT. |
| CMP-007 | Качество публикуется вместе с wall time, cell/element updates/s и memory footprint. |
| CMP-008 | Интерполяция между сетками консервативна либо её ошибка оценена отдельно. |

### 4.5. Требования к отчёту и развитию DG

| ID | Требование |
|---|---|
| RPT-001 | `docs/report.tex` является основной РПЗ НИРС; Markdown-версия не противоречит ей. |
| RPT-002 | В отчёте разделены «реализовано», «проверено», «измерено», «гипотеза» и «план». |
| RPT-003 | Заголовок, аннотация и выводы не утверждают гарантию `div B`, консервативность AMR, готовность GPU-порта или ускорение без пройденного gate. |
| RPT-004 | Каждая таблица/рисунок имеет case ID, mesh/resolution, physical time, solver profile и источник данных. |
| RPT-005 | Структура РПЗ: постановка; методы legacy/new; методика V&V; текущие результаты; ограничения; план AMR/GPU; выводы. |
| EXT-001 | Добавить ADR с границами physics, spatial operator, time integration, mesh/backend и diagnostics. |
| EXT-002 | FV остаётся отдельной реализацией spatial operator; DG не добавляется условными ветками в FV kernels. |
| EXT-003 | Общий state/flux API покрывается unit tests до рефакторинга под DG. |
| EXT-004 | Любой подготовительный рефакторинг сохраняет численный FV regression baseline. |
| EXT-005 | Реализация DG начинается только с отдельного ТЗ: basis/quadrature, numerical flux, limiter/positivity, divergence control, CFL и convergence suite. |

## 5. Нефункциональные требования

- Базовая acceptance-платформа — Linux x86-64; macOS является developer
  platform. Версии compiler/MPI/AMReX/CUDA зафиксированы.
- Ни один canonical test не содержит NaN/Inf или отрицательных `rho`, `p`.
  Floors не считаются доказательством устойчивости.
- Release benchmark включает минимум один warm-up и пять измеряемых повторов;
  сохраняются все времена, median и разброс.
- Compute-only runs отключают plot/diagnostic I/O; I/O измеряется отдельно.
- Strong scaling: `S_p=T_1/T_p`, `E_p=S_p/p`; weak scaling:
  `E_p=T_1/T_p` при постоянной нагрузке на ресурс.
- Для GPU отдельно измеряются kernel, transfers, MPI и end-to-end; используются
  Nsight Systems и Nsight Compute.
- Целевое ускорение утверждается после профилирования baseline. Непроверенное
  число ускорения не является контрактным требованием.
- Numerical change, optimization и dependency update — разные review units.

## 6. Начальные критерии приёмки

Пороги для double precision калибруются один раз на baseline и далее
версионируются. Ослабление требует численного обоснования.

| Область | Начальный gate |
|---|---|
| Kernel CPU | Обычно relative error `<= 1e-12`. |
| Legacy serial parity | Bitwise при неизменном порядке операций, иначе normalized `Linf <= 1e-12`. |
| Legacy OpenMP/CUDA parity | Short-run normalized `Linf <= 1e-10`; quality metrics меняются не более чем на 1%. |
| New CPU/GPU parity | Short-run normalized `Linf <= 1e-10`; end-time threshold калибруется по reductions. |
| Positivity | `rho>0`, `p>0`; canonical tests не используют скрытые floors. |
| New uniform CT | `dx*||div B||inf/max(||B||inf,Bref) <= 1e-12`. |
| Legacy `div B` | Absolute/scaled нормы опубликованы; лимит роста фиксируется после корректной инициализации baseline. |
| Uniform conservation | Periodic mass/energy relative defect `<= 1e-10`. |
| AMR conservation | После reflux relative defect `<= 1e-9` на regression case. |
| Alfvén convergence | Новый код: наблюдаемый порядок `p >= 1.8`; legacy first-order измеряется отчётно. |
| Brio–Wu | Case-specific thresholds по error/TV/front metrics после выбора reference branch. |

Для Brio–Wu «немонотонность» не сводится к одному графику: физическая
составная волна и неединственность решения отделяются от сеточного ringing.

## 7. Обязательная тестовая матрица

| Тест | Legacy | AMReX uniform | AMReX AMR | Назначение |
|---|---:|---:|---:|---|
| Constant state | да | да | да | preservation |
| HLLD states/degeneracies | да | да | — | kernel correctness |
| Brio–Wu | да | да | да | разрывы/немонотонность |
| Circular Alfvén wave | да | да | да | порядок/dissipation |
| Magnetic loop | да | да | да | CT/`div B` |
| Orszag–Tang | да | да | да | сложные волны |
| Rotor | желательно | да | да | magnetic braking |
| MHD blast | historical | да | да | robustness/positivity |
| Regrid conservation | — | — | да | coarse/fine correctness |
| MPI decomposition parity | — | да | да | decomposition independence |
| CPU/GPU parity | да | да | да | backend correctness |

## 8. Архитектура поставки

```text
MHD-DUAL
├── mhd2d-legacy repository
│   ├── solver + mesh adapter
│   └── vkr / corrected / opt_cpu / opt_cuda
├── mhd-amrex repository
│   ├── uniform + AMR solver
│   └── CPU / MPI / OpenMP / CUDA
└── shared validation contract
    ├── cases + mesh/reference manifests
    ├── adapters + converters + metrics
    └── raw -> summary -> figures -> report
```

Shared contract допустимо хранить в `mhd-amrex`, но без зависимости от
внутренних типов AMReX.

## 9. План работ

### P0 — baseline и управление требованиями

Зафиксировать source-of-truth legacy, immutable tag, dirty patch и mesh
generator; зафиксировать commits нового кода/AMReX; утвердить canonical test
cards, Brio reference, hardware matrix и requirements-to-test traceability.

### R0 — исправление отчёта для научного руководителя

Провести claim audit `REPORT.md`/`report.tex`; убрать или маркировать
неподтверждённые гарантии; описать legacy и новую схему, текущий evidence и
ограничения; добавить программу экспериментов и гипотезу H1; привести документ
к сокращённой структуре РПЗ НИРС. Текущие диагностические числа пометить как
не-benchmark. Выход: компилируемый draft PDF и таблица claim-to-evidence.

### L0 — восстановление legacy

Очистить сборку без изменения математики, восстановить Netgen pipeline,
разделить source/input/output, добавить config/physical time и воспроизвести
Brio, Alfvén и magnetic loop из ВКР. Выход: `legacy_vkr` и reproduction report.

### L1 — корректность legacy

Добавить unit/regression tests HLLD, conversions, geometry и EMF; исправить
CPU/CUDA drift; исправить divergence/integral diagnostics; добавить invalid-state
dumps. Выход: `legacy_corrected` и bug-fix ledger.

### L2 — CPU/OpenMP legacy

Снять serial profile/roofline estimate, исправить OpenMP, убрать повторные
allocation/copies, оптимизировать layout только по данным профиля и провести
serial/OpenMP scaling с parity. Выход: `legacy_opt_cpu`.

### L3 — CUDA legacy

Сначала воспроизвести partial-offload ВКР, затем обеспечить HLLD parity,
persistent allocations, resident hot path и минимальные host round trips;
провести Nsight-анализ. Выход: `legacy_opt_cuda` и CPU/GPU report.

### N0 — инженерный baseline AMReX

Стабилизировать build matrix, добавить CTest, manifests, timing regions и
контроль output schedule. Выход: воспроизводимый CPU release.

В рамках N0 также создать минимальный ADR модульных границ для будущего DG.
Не выполнять большой рефакторинг до появления regression baseline.

### N1 — uniform-grid корректность

Добавить HLLD fallback/counters, проверить limiter и SSP-RK2 convergence,
согласовать flux/CT EMF и выполнить canonical suite.

### N2 — Brio–Wu

Утвердить reference branch, реализовать versioned front identification,
посчитать error/TV/extrema/front metrics и выполнить ablation/convergence.

### N3 — консервативный AMR

Реализовать gas reflux, restriction/prolongation/sync staggered `B`, добавить
regrid/conservation tests и multilevel canonical runs.

### N4 — CPU/MPI/OpenMP performance

Профилировать compute/communication/regrid/diagnostics/I/O, устранить serial
bottlenecks и провести single-/multi-node strong/weak scaling с binding data.

### N5 — CUDA AMReX

Перевести hot path и служебные loops на performance-portable kernels, устранить
host-only captures/неявные sync, пройти общий regression suite и снять Nsight
profiles/GPU scaling.

### C0/C1/C2 — общий эксперимент

Создать schema/adapters/converters и эквивалентную меру resolution; выполнить
equal-resolution, equal-cost, equal-error comparison; затем кластерную campaign,
проверку raw data и автоматическую генерацию отчёта.

## 10. Очерёдность и контрольные ворота

```text
P0
├─ R0 draft -------------------------------> C2 final report
├─ L0 -> L1 -> L2 -> L3 ─┐
├─ N0 + DG ADR -> N1 -> N2 -> N3 -> N4 -> N5
└─ C0 --------------------┴-> C1 -> C2
```

| Gate | Условие |
|---|---|
| G-P0 | Утверждены commits, test cards, reference и hardware plan. |
| G-R0 | РПЗ собирается; все claims классифицированы, неподтверждённые гарантии удалены. |
| G-L0 | `legacy_vkr` clean-build и минимум три воспроизведённых теста ВКР. |
| G-L1 | Corrected suite проходит; numerical deltas перечислены. |
| G-L2 | CPU optimization проходит parity и даёт статистически значимый эффект. |
| G-L3 | CUDA parity; опубликованы kernel и end-to-end timings. |
| G-N0 | Clean build + CTest + manifests без ручных шагов. |
| G-N1 | Uniform numerical gates выполнены. |
| G-N2 | Brio metrics и ablation выполнены на утверждённом reference. |
| G-N3 | AMR conservation/regrid gates выполнены. |
| G-N4 | CPU scaling campaign воспроизводима. |
| G-N5 | CUDA regression/profiling gates выполнены. |
| G-C1 | Завершены три режима прямого сравнения. |
| G-C2 | Cluster raw data проверены; отчёт генерируется автоматически. |
| G-EXT0 | ADR для FV/DG-расширения принят; текущий FV regression не изменён. |

## 11. Трудоёмкость и команда

| Направление | Оценка, чел.-нед. |
|---|---:|
| Baseline, требования, shared schema | 1–2 |
| Legacy recovery/reproduction | 2–4 |
| Legacy correctness/tests | 2–3 |
| Legacy CPU/CUDA optimization | 3–6 |
| AMReX uniform correctness/Brio | 3–5 |
| AMReX conservative AMR | 3–5 |
| AMReX CPU/GPU performance | 4–8 |
| Cross-validation/cluster/report | 3–5 |
| **Последовательное effort** | **18–38** |

Для команды из трёх независимых потоков ориентир — 9–16 календарных недель
после получения внешних данных и кластера. Это planning range, не обещанная
дата. Роли: technical lead, numerical V&V, legacy, AMReX/AMR,
GPU/performance, cluster/reproducibility и report owner.

AI-агент получает один work package с входным commit, profile, test card и
acceptance gate. Один интегратор согласует формулы/config. Разные агенты не
редактируют одновременно одни numerical kernels.

## 12. Риски

| Риск | Управление |
|---|---|
| Неоднозначный legacy snapshot | Immutable tag + dirty patch + archive checksum |
| Bug fix смешан с optimization | Разные profiles/commits + parity suite |
| Неединственность Brio–Wu | Выбрать regular/compound reference до финального gate |
| Старый GPU путь не активен | Staged reproduction; не выдумывать GPU-числа |
| Нет AMR reflux | G-N3 блокирует claims о global conservation |
| CPU loops в CUDA build | Device coverage audit + profiler evidence |
| Несопоставимые сетки/стоимость | Три режима сравнения + явная мера `h` |
| Подтверждение желаемой гипотезы подменяет исследование | H1 имеет falsification criteria; публикуются и отрицательные результаты |
| Преждевременный DG-рефакторинг ломает FV | Сначала tests/ADR; DG — отдельное ТЗ и change stream |
| Архив ~60 GB без provenance | Только forensic data; новые runs с manifests |
| Нет кластера | Заранее job scripts; G-C2 остаётся открытым |

## 13. Поставка и Definition of Done

Release bundle содержит два tagged source tree, build/toolchain configs,
canonical cases/mesh manifests, CTest/regression/parity suites, benchmark и
cluster scripts, `benchmarks/raw`, `summary`, `figures`, Nsight summaries,
traceability matrix, reproduction/verification/performance/final reports.

Программа завершена, когда:

1. оба продукта собираются и запускаются из clean checkout;
2. legacy имеет раздельные `vkr/corrected/opt_cpu/opt_cuda` профили;
3. результаты ВКР воспроизведены в утверждённом допуске либо расхождение строго
   объяснено;
4. новый solver проходит uniform, AMR, MPI и CPU/GPU gates;
5. AMR-консервация доказана с reflux и контролем `div B`;
6. Brio–Wu исследован машинными метриками на выбранном reference;
7. завершены equal-resolution/cost/error comparisons;
8. scaling содержит raw runs, статистику и environment data;
9. каждый вывод отчёта трассируется до конфигурации и данных;
10. нет заявлений «реализовано/ускорено/точнее» без прошедшего gate;
11. принят ADR расширения FV/DG, не изменивший текущий FV regression baseline.

## 14. Внешние решения до финальной приёмки

- выбрать официальный legacy source: локальный `9d0f60e`, архивный `3023cb`
  либо другой commit;
- выбрать regular или compound reference Brio–Wu;
- предоставить точную работу Бисикало–Жилкина;
- предоставить cluster scheduler/account/partition/modules и hardware;
- определить GPU architectures и допустимые CUDA/MPI versions;
- определить права публикации исходников и результатов ВКР.
- предоставить шаблон/методические требования кафедры к РПЗ НИРС, если они
  обязательны; до этого используется сокращённая инженерно-научная структура.

До получения этих данных можно готовить код, тесты и job scripts, но нельзя
заполнять итоговые comparative/scaling tables предполагаемыми числами.
