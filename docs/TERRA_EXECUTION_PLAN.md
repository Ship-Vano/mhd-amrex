# План исполнения ТЗ агентом Terra

Статус: рабочая декомпозиция на фазы, 2026-08-26.

Этот документ превращает [`PROGRAM_PLAN.md`](PROGRAM_PLAN.md) и
[`TECHNICAL_SPECIFICATION.md`](TECHNICAL_SPECIFICATION.md) в последовательность
ограниченных заданий для GPT-5.6 Terra Ultra. Он не заменяет математические
требования этих документов. Если краткая карточка фазы расходится с детальным
ТЗ, действует более строгий критерий.

## 1. Как использовать план

Один запуск/контекст Terra должен выполнять одну фазу или один явно выделенный
подэтап. Не выдавать агенту запрос «заверши весь проект»: он смешает baseline,
численную коррекцию, оптимизацию и оформление результатов.

Для каждой фазы владелец передаёт:

- входной Git commit и ожидаемую ветку;
- разрешённые каталоги/файлы;
- доступные внешние данные и оборудование;
- идентификатор фазы и её gate;
- запрет на переход к следующей фазе без проверки.

Terra обязан завершить фазу одним из состояний:

- `PASS` — gate доказан командами и артефактами;
- `PARTIAL` — полезная часть готова, но gate не пройден;
- `BLOCKED` — отсутствует внешнее решение/доступ;
- `FAILED` — проверка выявила дефект, требующий нового work package.

Статус `PASS` нельзя выводить только из компиляции или визуального графика.

## 2. Общий контракт каждой фазы

Перед изменениями агент обязан:

1. прочитать `AGENTS.md` и относящиеся разделы основного ТЗ;
2. выполнить `git status`, `git log -1`, записать входной commit/dirty-state;
3. воспроизвести релевантную сборку и существующие тесты;
4. сформулировать гипотезу изменения и измеримый acceptance criterion;
5. перечислить файлы, которые предполагается изменить.

Во время работы:

- один commit решает одну проблему;
- test/diagnostic должен появляться до или вместе с исправлением;
- numerical fix не смешивается с optimization;
- dependency upgrade не смешивается с numerical change;
- пользовательские plotfile, архивы и внешние legacy-каталоги не удаляются;
- любые новые числа сохраняются как raw artifact с manifest;
- отрицательный результат не скрывается.

В handoff агент обязан указать:

```text
Phase:
Input commit:
Output commit(s):
Status: PASS | PARTIAL | BLOCKED | FAILED
Changed:
Verified by:
Raw artifacts:
Measured facts:
Inferences:
Open risks/blockers:
Allowed next phase:
```

## 3. Карта фаз

| Фаза | Назначение | Основная зависимость | Параллельность |
|---|---|---|---|
| T00 | Baseline и decision register | текущий репозиторий | первая |
| T01 | Исправленная РПЗ НИРС, draft | T00 | параллельно T02/T03 |
| T02 | Build/CTest/manifest foundation | T00 | параллельно T01/T03 |
| T03 | Legacy provenance и `legacy_vkr` | T00 + решение владельца | отдельный source tree |
| T04 | Uniform numerical correctness | T02 | critical path |
| T05 | Legacy corrected/parity baseline | T03 | параллельно T04 |
| T06 | Brio–Wu metrics и диагноз | T04 + reference; T05 для сравнения |
| T07 | Conservative AMR | T04 | не параллельно с изменениями CT/HLLD |
| T08 | Cross-solver quality comparison | T05 + T06 + T07 |
| T09 | Timing и CPU benchmark foundation | T02 + correctness gates |
| T10 | Legacy CPU/CUDA optimization | T05 + T09 | отдельный legacy tree |
| T11 | AMReX CPU/MPI/OpenMP scaling | T07 + T09 |
| T12 | AMReX CUDA port и parity | T07 + T09 |
| T13 | Cluster campaign | T10/T11/T12 + cluster access |
| T14 | Финальная РПЗ и release evidence | T08 + доступные performance gates |
| TDG0 | ADR и seam для будущего DG | T04 regression baseline | без реализации DG |

`T01` создаёт честный промежуточный документ сейчас; `T14` заменяет плановые
разделы измеренными результатами после завершения вычислительных фаз.

## 4. Critical path

```text
T00 -> T02 -> T04 -> T06 -> T07 -> T08 -> T14
  │      │                           ↑
  │      └-> T09 -> T11/T12 -> T13 -┘
  ├-> T01 --------------------------> T14
  └-> T03 -> T05 -> T10 -> T13 -----┘
                 └-----------> T08

T04 -> TDG0
```

До T04 запрещены сильные утверждения о корректности uniform solver. До T07
запрещены утверждения о глобальной AMR-консервативности. До T08 нельзя заявлять
превосходство AMReX над legacy. До T12 нельзя называть сборку GPU-портом.

## 5. Карточки фаз

### T00 — baseline и решения

**Цель:** создать однозначную точку отсчёта и список внешних решений.

**Работы:**

- проверить clean/dirty state, build environment и существующие быстрые тесты;
- сверить `PROGRAM_PLAN`, детальное ТЗ, legacy/external manifests;
- создать/обновить requirements-to-test и decision register;
- перечислить неподтверждённые claims текущего отчёта;
- зафиксировать source candidates legacy без выбора за владельца;
- зарегистрировать blockers: Brio reference, статья, cluster, НИРС template.

**Артефакты:** baseline manifest, decision register, status table.

**Gate T00:** ни один входной commit, внешний dataset или открытое решение не
остались неявными.

**Запрещено:** менять численную схему, оптимизировать, заполнять отсутствующие
результаты предполагаемыми числами.

### T01 — draft РПЗ для научного руководителя

**Цель:** привести `report.tex` и `REPORT.md` к честному текущему состоянию.

**Работы:**

- создать claim-to-evidence matrix;
- убрать неподтверждённую «гарантию» из названия и выводов;
- описать реализованные HLLD/MUSCL/CT/SSP-RK2 и явно отделить их проверку;
- обозначить отсутствие gas reflux, verification suite, реального GPU-порта и
  воспроизводимого scaling;
- добавить H1, методику её будущей проверки и roadmap;
- пометить единичные времена как diagnostic, не benchmark;
- собрать PDF и проверить ссылки/рисунки/нумерацию.

**Gate T01:** PDF собирается; каждое существенное утверждение имеет evidence
status; документ можно передать руководителю без завышения зрелости.

**Запрещено:** генерировать фиктивные таблицы, объявлять план результатом,
параллельно рефакторить solver ради текста отчёта.

### T02 — инженерная основа тестов и запусков

**Цель:** сделать каждый следующий результат воспроизводимым.

**Работы:**

- CMake presets и CTest registration;
- run manifest с commit/config/dependencies/hardware;
- строгая JSON validation;
- benchmark mode без некритического I/O;
- структура `tests/`, `benchmarks/raw|summary|figures`, scripts;
- устранение двойного финального plotfile отдельным малым change set.

**Gate T02:** clean configure/build/CTest; намеренная поломка тестируемого
kernel обнаруживается; повторный запуск сохраняет эквивалентный manifest.

**Запрещено:** менять HLLD/MUSCL/CT, публиковать новые performance claims.

### T03 — provenance и воспроизведение legacy

**Цель:** получить неизменяемый профиль `legacy_vkr`.

**Внешний вход:** подтверждение владельца, какой commit является официальным;
привязка mesh/config/results либо решение перезапустить задачи.

**Работы:**

- immutable tag и archive checksum;
- документированная clean build без абсолютных путей;
- воспроизводимый Netgen mesh pipeline;
- config/physical time/run manifest без изменения математики;
- повтор Brio–Wu, Alfvén и magnetic loop;
- comparison с опубликованными данными ВКР с явными допусками.

**Gate T03:** три canonical legacy cases воспроизводятся; источник каждого
файла и числа доказан.

**Запрещено:** исправлять формулы в `legacy_vkr`, считать файлы на Elements
доказательством точного запуска без manifest.

### T04 — uniform-grid numerical correctness

**Цель:** доказать корректность нового solver без AMR.

**Подэтапы, каждый отдельным commit:**

1. state/conversion/wave-speed unit tests;
2. HLLD consistency/branches/degeneracies;
3. validity predicate, HLLE fallback и counters;
4. limiter tests и reconstruction invariants;
5. CT constant-state/1D/divergence tests;
6. SSP-RK2 Butcher/convergence tests;
7. Alfvén spatial-temporal convergence.

**Gate T04:** все CPU Release CTest проходят; стандартные cases finite и
positive; fallback/floor policy наблюдаема; Alfvén достигает принятого порядка;
uniform CT выполняет заданную норму `div B`.

**Запрещено:** добавлять AMR reflux, GPU-порт или performance optimization в те
же commits; скрывать отрицательные состояния floor без счётчика.

### T05 — `legacy_corrected` и parity baseline

**Цель:** исправить известные legacy-дефекты, не уничтожив исторический профиль.

**Работы:**

- tests для conversions/HLLD/geometry/EMF;
- исправление CPU/CUDA HLLD drift только в corrected profile;
- `max(abs(div B))` и scaled norm;
- positivity/finiteness diagnostics и failure dump;
- numerical-delta ledger между `vkr` и `corrected`.

**Gate T05:** `legacy_vkr` не изменён; corrected tests проходят; каждый delta
объяснён и имеет regression.

**Запрещено:** оптимизировать data layout или OpenMP/CUDA в этом change set.

### T06 — Brio–Wu: метрики и причина немонотонности

**Цель:** заменить визуальное суждение машинно проверяемым диагнозом.

**Внешний вход:** выбранная regular либо compound/non-regular reference branch.

**Работы:**

- зафиксировать одинаковые physics/BC/CFL/times;
- reference projection и versioned front annotations;
- `L1/L2/Linf`, excess TV, overshoot/undershoot, extrema count;
- положение/толщина каждого фронта;
- grid convergence и ablation first-order/MUSCL/SSP-RK2/full CT;
- legacy comparison после T05;
- поставить диагноз сначала на uniform grid; AMR-вариант исследовать только
  после T07 и не объявлять адаптацию лечением заранее.

**Gate T06:** каждый экстремум классифицирован как reference structure,
discretization artifact или неразрешённая неопределённость; метрики
воспроизводятся одной командой.

**Запрещено:** выбирать reference после просмотра желаемого результата,
оценивать монотонность по одному графику.

### T07 — conservative AMR

**Цель:** закрыть correctness-критичный coarse/fine defect.

**Подэтапы:**

1. conservation diagnostic, который демонстрирует текущий дефект;
2. gas flux register/reflux для mass/momentum/energy;
3. face-`B` restriction/prolongation и CT synchronization audit;
4. regrid positivity/`div B` checks;
5. periodic/closed regression и uniform-fine comparison;
6. MPI decomposition check.

**Gate T07:** AMR conservation и normalized `div B` удовлетворяют заранее
заданным tolerances; test падает при отключённом reflux; error AMR не хуже
uniform coarse и приближается к uniform fine.

**Запрещено:** одновременно менять HLLD/limiter, заявлять conservation только
по массе либо по одному уровню.

### T08 — прямое сравнение качества

**Цель:** проверить accuracy-часть H1 без performance cherry-picking.

**Работы:**

- canonical adapters legacy/new;
- эквивалентная мера `h` и conservative comparison projection;
- equal resolution, equal cost, equal error;
- uniform AMReX, fixed AMR и dynamic AMR отдельно;
- error/front/divergence/conservation/memory таблицы;
- Pareto frontier и объяснение результатов по задачам.

**Gate T08:** все таблицы строятся из raw/summary; ни одно улучшение не
приписано AMR без ablation; H1 получает status supported/partial/refuted по
заранее заданным cases.

### T09 — timing и benchmark foundation

**Цель:** обеспечить честное измерение до оптимизаций.

**Работы:** regions `init/compute/communication/regrid/diagnostics/I/O/E2E`,
warm-up + минимум пять repeats, median/MAD, MPI max wall time, affinity/binding,
cell updates/s и machine manifest. Compute-only и I/O experiments разделены.

**Gate T09:** повторяемость на локальной машине находится в принятом диапазоне;
raw records полны; compile/FetchContent/I/O не попадают в compute time.

### T10 — оптимизация legacy

**Цель:** получить `legacy_opt_cpu` и `legacy_opt_cuda` без изменения качества.

**Порядок:** profile → CPU/OpenMP → parity → CUDA historical reproduction →
persistent device data/hot path → CPU/GPU parity → Nsight.

**Gate T10:** ускорение статистически подтверждено относительно corrected
baseline; quality delta опубликован; kernel и end-to-end времена разделены.

**Запрещено:** переносить математические исправления в perf commit, сравнивать
один CPU core с целым GPU как основной результат.

### T11 — AMReX CPU/MPI/OpenMP scaling

**Цель:** проверить CPU performance-часть H1.

**Работы:** serial/OMP/MPI/hybrid single-node; затем strong/weak multi-node;
fixed AMR и dynamic AMR отдельно; rank/thread binding и load balance metrics.

**Gate T11:** опубликованы raw repeats, `S_p`, `E_p`, cells per rank и breakdown;
точность расчёта не изменилась относительно correctness baseline.

### T12 — AMReX CUDA port

**Цель:** получить реальное device execution, а не CUDA-compatible headers.

**Подэтапы:** build/architecture; hot-path `ParallelFor/ParReduce`; init/BC;
`ComputeDt/MaxDivB`; tagging/regrid; derived/output; CPU/GPU regression;
compute-sanitizer; Nsight Systems/Compute; multi-GPU MPI.

**Gate T12:** timed hot path не содержит `LoopOnCpu`; CPU/GPU parity и общий
regression suite проходят; kernel/transfer/MPI/E2E измерены.

**Запрещено:** host-only capture/device dereference, performance claim по
асинхронному времени без корректной синхронизации.

### T13 — cluster campaign

**Цель:** получить финальные CPU/GPU scaling data.

**Внешний вход:** scheduler, account/partition, limits, modules, node topology,
GPU-aware MPI и profiling permissions.

**Работы:** reproducible build scripts, job arrays с независимыми run IDs,
strong/weak CPU, GPU scaling, profile representative cases, raw validation.

**Gate T13:** каждый summary row восстанавливается из raw runs и environment;
failed/outlier runs сохранены и классифицированы, а не удалены молча.

### T14 — финальная РПЗ и release evidence

**Цель:** заменить roadmap draft T01 доказанными результатами.

**Работы:** обновить claim matrix; сгенерировать таблицы/рисунки из summary;
изложить H1 per case; разделить результаты и ограничения; добавить AMR/GPU/
scaling только для пройденных gates; выпустить reproduction instructions.

**Gate T14:** PDF собирается одной документированной командой; цифры совпадают
с summary; conclusions трассируются к raw data; руководитель видит как
положительные, так и отрицательные результаты.

### TDG0 — архитектурный ADR для будущего DG

**Цель:** оставить расширяемую основу после доказательства FV baseline.

**Работы:** dependency map; ADR границ `physics/spatial/time/backend/diagnostics`;
минимальный seam только при наличии текущего потребителя; unit contracts общего
state/physical-flux API; отдельный будущий backlog DG basis/quadrature/surface
flux/limiter/positivity/divergence control/CFL/convergence.

**Gate TDG0:** ADR принят; все FV regression tests дают прежний результат;
нет DG-заглушек, runtime `if (dg)` в FV kernels и speculative framework.

## 6. Разрешённая параллельность Terra Ultra

После T00 допустимы максимум три независимых потока:

1. **Report/V&V data:** T01 и подготовка reference literature/data без изменения
   numerical kernels.
2. **AMReX correctness:** T02 → T04 → T06/T07; один владелец HLLD/CT/AMR.
3. **Legacy:** T03 → T05 → T10 в отдельном source tree.

T09 может выполняться параллельно T06 только после стабилизации интерфейса
diagnostics. T11 и T12 разделяются лишь после T07 и не меняют общие numerical
kernels без возврата к correctness owner. T08 и T14 выполняет интегратор.

Нельзя параллельно:

- исправлять один HLLD/CT kernel в двух ветках;
- менять case normalization и одновременно собирать comparative data;
- оптимизировать и менять acceptance tolerance;
- писать итоговые выводы до freeze summary;
- выполнять DG-рефакторинг одновременно с построением FV baseline.

## 7. Шаблон запроса для Terra

```text
Выполни только фазу <Txx> из docs/TERRA_EXECUTION_PLAN.md.

Input commit: <sha>
Разрешённый scope: <paths/subtask>
Доступные данные/оборудование: <list>
Запрещено: переходить к следующей фазе; смешивать numerical/performance changes;
изменять tolerances без обоснования; выдумывать результаты.

Сначала прочитай AGENTS.md, PROGRAM_PLAN.md и относящиеся разделы
TECHNICAL_SPECIFICATION.md. Проверь baseline. Сформулируй гипотезу и gate.
Реализуй минимальный связный change set, добавь regression, выполни проверки.

Заверши handoff по шаблону из раздела 2. Если внешний blocker не мешает
подготовительной работе, выполни её и поставь PARTIAL. Commit создавай только
после diff review и успешных релевантных проверок.
```

## 8. Начальная очередь

Рекомендуемый порядок ближайших запусков Terra:

1. `T00` — короткая актуализация baseline/decisions.
2. `T01` — исправление РПЗ для научного руководителя.
3. `T02` — CTest/manifests/benchmark foundation.
4. `T03` — после решения владельца о legacy source; до решения подготовить
   только provenance plan.
5. `T04` — последовательные небольшие kernel/correctness change sets.
6. `TDG0` — только после стабильного T04 regression baseline.

На текущем состоянии нельзя начинать с T07, T10, T11 или T12: без T02/T04
полученное ускорение не будет численно доказанным.
