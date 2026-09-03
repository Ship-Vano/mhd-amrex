# Состояние legacy_corrected: срез T05

Дата среза: 2026-08-31. Статус: **PARTIAL**. Этот файл создан как
evidence-bound handoff при досрочной фиксации работы: он сохраняет полезные
результаты, но не переименовывает незакрытую проверку в завершённую фазу.

## Что зафиксировано

Неизменяемый исторический source tree
`/Users/ivansamanov/Documents/MHD2D` проверен чистым на commit
`9d0f60ea8576fac5d6f28c4dec142236d76131d6`. Он не модифицировался. Новый
CPU-профиль живёт только как versioned overlay
[`legacy/patches/0001-legacy-corrected-physics.patch`](../legacy/patches/0001-legacy-corrected-physics.patch),
SHA-256 `7d7314c69c0813a989a6bc13591d38f5bce4966d69c3be7a1c8317818477e6bf`.

На чистом clone указанного commit были выполнены `git apply --check
--whitespace=error`, Release configure/build GCC 15.2 с OpenMP 4.5 и
`ctest --output-on-failure`: `1/1` тест пройден. В CTest входят conversion
round-trip, HLLD/HLLE, Brio--Wu flux, 1D CP-Alfvén, корректность геометрии и
ghost mapping, RT0 на действительно неправильных треугольниках, ориентация CT
для CW-входа, discrete curl/divergence и CFL для sliver-треугольника.

## Сохранённые измерения

Сырые VTU, mesh, config, log и manifest намеренно игнорируются Git и сохранены
в `benchmarks/raw/legacy_corrected/`. Компактная машиночитаемая сводка лежит в
[`benchmarks/summary/legacy_corrected_t05_state.json`](../benchmarks/summary/legacy_corrected_t05_state.json).
Оба завершённых run получили `quality_gate.status = pass`: процесс завершился
нормально, достигнут заданный физический момент, итоговое состояние конечно и
положительно; для периодических задач прошёл balance полной энергии.

| Case | Неструктурированная сетка | Измеренный итог | Физическая интерпретация |
|---|---:|---|---|
| CP-Alfvén, Tóth 30°, `gamma=5/3`, `CFL=0.1`, `t=1` | 3348 треугольников, requested `maxh=0.04`, min angle `37.74°` | `rho_min=0.904806`, `p_min=0.0942713`, fallback `0`, `max|div B|=6.089e-13`, scaled `3.379e-15`, residual полной энергии `5.106e-15`; средняя relative-L1 wave metric `0.301561` | расчёт физически допустим и дискретно соленоидален в измеренной норме, но ошибка велика. Это не совпадение с таблицей второго порядка Tóth. |
| Field loop, Athena-compatible inclined geometry, `gamma=5/3`, `CFL=0.1`, один период `t=8/sqrt(3)` | 828 треугольников, requested `maxh=0.08`, min angle `38.97°` | `rho_min=1`, `p_min=1.000000224`, fallback `0`, `max|div B|=1.400e-16`, scaled `1.661e-15`, residual энергии `-1.155e-14`; `E_B(t)/E_B(0)=0.0473140`, return-B relative L1 `1.23905` | состояние устойчиво, но перенос магнитной петли на этой грубой сетке очень диффузен. Результат отрицательный по качеству и сохранён именно поэтому. |

Для CP отдельный 1D HLLD regression без неструктурированной CT-связки даёт
relative L1 `0.425793`, `0.242329`, `0.129581` при `N=32,64,128`; ошибка
монотонно уменьшается, fallback не срабатывает. Это проверяет согласованность
исправленного потока, но не является доказательством второго порядка: legacy
остается first-order Euler/FV методом.

Завершённые raw runs использовали предыдущий overlay
`32f6260be09a3db87091be23c98f48fba218ee2ef9ae70df6b92dc793b29bc5a`.
После них в текущий overlay добавлен только conversion round-trip regression;
`git diff` между двумя overlays содержит conversion regression и механическое
удаление trailing whitespace; production solver/geometry semantics не менялись.
Поэтому числа сохранены как валидные измерения именно того run,
но для финальной post-commit кампании их необходимо повторить с SHA текущего
overlay, а не выдать за результаты новой ревизии.

Два начатых после этого запуска — `magnetic_loop_athena_h004_t05final` и
`brio_wu_h0025_t05final` — были прерваны по запросу владельца во время solver
stage. Их неполные ignored-директории оставлены на диске; ни одна величина из
них не входит в таблицу выше.

## Numerical-delta ledger

| Delta corrected CPU-профиля | Причина и физический смысл | Regression/наблюдение | Граница утверждения |
|---|---|---|---|
| `HLLD_flux_corrected` + наблюдаемый HLLE fallback | промежуточные HLLD состояния и сигнал-скорости проверяются на конечность/положительность; недопустимый промежуточный HLLD не должен распространять NaN | `F(U,U)=F_phys`, обе supersonic ветви, Brio--Wu, `Bn=0`, near-vacuum fallback и rejection invalid state | это CPU flux safeguard, а не глобальная гарантия positivity |
| Primitive/conservative conversion regression | давление должно вычитать кинетическую и магнитную энергии ровно один раз | round-trip `rho,u,v,w,p` в CTest | покрыт core conversion, не все исторические I/O paths |
| CFL для произвольного треугольника | `h_min` не контролирует спектральный объём sliver-клетки | CTest сверяет `CFL*A/sum(lambda*L)` и проверяет отсутствие скрытого нижнего floor | применено к corrected CPU path; это не оптимизация |
| Граница--ghost | один boundary element может иметь два граничных ребра; отображение только element→ghost теряло угол | взаимно-однозначные edge↔ghost maps на irregular mesh | corrected path требует text mesh |
| `div B` и CT orientation | signed maximum мог скрывать отрицательную ошибку; знак закона Фарадея зависит от edge tangent/normal orientation | `maxAbs`, flux и scaled norms; potential curl и CT update сохраняют divergence; CW text mesh канонизируется, отрицательная helper-ветвь тестируется | machine-level invariant для проверенной topology, не theorem для всех исторических binary meshes |
| RT0 reconstruction | cell `B` должен восстанавливаться из face-normal `B` со знаком относительно элемента, а не из порядка узлов | constant physical `B=(0.73,-0.41)` восстанавливается на неправильных треугольниках до `2e-12` | RT0 projection создаёт представление, не точное pointwise CP поле |
| Энергия после CT | подгонка `E` после RT0-реконструкции `B` искусственно меняла полную энергию | исправленный путь сохраняет консервативную `E`; отдельно логируется изменение магнитной энергии от reconstruction | CT/EMF остаётся low-order nodal averaging, не upwind CT |
| Failure diagnostics и manifests | exit code сам по себе не отличает физический failure | `physical_failure.json`, min `rho/p`, fallback count, divergence, balances, mesh/config/source hashes и quality gate в runner | это CPU run provenance, не benchmark |
| Карты 4 и 5 читают конфиг (CFL, `finalTime`) | `task_type` 4 (цилиндр) и 5 (ОТ) жёстко задавали `cflNum=0.5` и `finalTime`, поэтому манифест кейса игнорировался и прогон нельзя было воспроизвести по конфигу | обе карты управляются `solverConfig.json` так же, как 1/8/9; прогон `benchmarks/raw/legacy_corrected/rotor_128/` с quality gate | отрицательный/отсутствующий `cfl` сохраняет исторические умолчания |
| Бездивергентная инициализация карт 4 и 5 | `B_n` бралось сэмплированием в середине ребра; для непостоянного поля контурная сумма по треугольнику не равна нулю на произвольной сетке. Теперь `B_n` --- разность узловых значений `A_z`, а клеточное `B` --- RT0-проекция тех же граней | начальная невязка магнитного потока ОТ `1.0e-17`, цилиндра `1.0e-16`; `rho_min` цилиндра сместился `0.764 -> 0.743` | на структурной прямоугольно-треугольной сетке старая схема давала нуль случайно (слагаемые телескопируются); дельта видна на нерегулярной сетке |
| Состояния на рёбрах карт 4 и 5 | `initEdgeUs` только резервировался и оставался нулевым (нулевая плотность), хотя `runSolver` копирует его в `edgeUs` | обе карты заполняют `initEdgeUs` аналитическим состоянием, периодические карты --- через образ внутри домена | наблюдаемого влияния на итог не выявлено; устранена скрытая зависимость от порядка перезаписи |
| Отказ вместо тихой подмены на неисправленных путях | `gpu:true` выполнял пустой блок: решатель не запускался, но записывался неэволюционировавший VTU с кодом 0; `cylindrical:true` исполнял неисправленные ядра | оба пути возвращают код 2 с явным сообщением | это отказ от claim, а не перенос исправлений на GPU/осесимметрию |
| Диагностика физического отказа переживает исключение | `physical_failure.json` оставался пустым: исключение не перехватывалось, `std::terminate` не разматывал стек и деструктор `ofstream` не выполнялся | поток явно сбрасывается перед `throw`; `main` перехватывает и возвращает 3; в дамп добавлены центроид, `rho`, `p` | дамп по-прежнему фиксирует одну ячейку --- ту, что обнаружена первой |

## Научные ограничения и результаты аудита

- CP поставлен по 30-градусной геометрии Tóth; Athena CP page использует другой
  угол. Кроме того, таблица Tóth получена второпорядковой base scheme. Поэтому
  `0.301561` нельзя ставить рядом с её числами как прямое quality comparison.
  Корректные источники постановки: [Tóth (2000)](https://public.websites.umich.edu/~gtoth/Papers/Toth2000_divb.pdf)
  и [Athena CP-Alfvén test](https://www.astro.princeton.edu/~jstone/Athena/tests/cp-alfven-wave/cp-alfven.html).
- Статья Avdeeva--Lukin использует Raviart--Thomas/edge-normal поле, HLLD и
  треугольную staggered-сетку, но её тестовые области указаны как сетки из
  равносторонних треугольников. Нынешняя Netgen-сетка намеренно неправильная,
  поэтому нет claim о точном воспроизведении её картинок или чисел. См.
  [DOI статьи](https://doi.org/10.1088/1742-6596/1336/1/012026).
- `legacy_vkr` по-прежнему не воспроизвёл исторические результаты: T03 остался
  `FAILED`; archived VTU на `Elements` не имеет доказанной связки
  commit/config/run. Исторические числа ВКР не переписывались.
- CUDA source не менялся и не проверялся; runner явно использует `gpu=false`.
  Следовательно, нет CPU/CUDA parity, GPU-correctness или performance claim.
- Нет выбранной regular/compound Brio--Wu reference branch. Нет front/TV/
  overshoot classification; это строго следующий T06, а не результат T05.
- Не сделаны performance/scaling measurements. GCC/OpenMP здесь лишь build
  dependency, а не измеренное ускорение.
- Binary `World` намеренно отклонён corrected runner: для него отсутствуют
  новые edge↔ghost maps. Это явная ограниченная поддержка, не молчаливая
  несовместимость.

## Gate и безопасная точка продолжения

Техническая часть CPU T05 сохранена: immutable `legacy_vkr`, применимый
overlay, passing CTest и delta ledger есть. Полный T05 handoff остаётся
**PARTIAL**, поскольку нет новых canonical raw runs для SHA текущего overlay,
CUDA parity не проверена, а Brio--Wu физически не классифицирован. Это не
мешает другому агенту начать с этого commit, но он должен сначала повторить
CP/field-loop/Brio с текущим SHA и сохранить manifests; затем отдельным
решением переходить к T06 только после выбора reference branch.

Воспроизводимая команда для нового CPU-run:

```sh
python3 -B scripts/run_legacy_corrected.py \
  --source /Users/ivansamanov/Documents/MHD2D \
  --case cp_alfven \
  --artifact-dir benchmarks/raw/legacy_corrected/<new-run-id> \
  --compiler /opt/homebrew/bin/g++-15 --omp-threads 1
```
