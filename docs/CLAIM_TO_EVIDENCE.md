# Claim-to-evidence matrix for the NIRS report

Version: 2026-08-31 (updated after T04 / T05-CPU / T08-partial / literature
validation). Status values have a strict meaning:

- `implemented` — code is present; not a numerical proof;
- `verified` — asserted by an automated test (CTest / runner quality gate);
- `measured` — a diagnostic observation, not a benchmark or full verification;
- `provisional` — measured, but against a reference that is not yet independent;
- `hypothesis` — research hypothesis, criterion defined, not settled;
- `not established` — needs a future phase / external decision;
- `blocked` — needs external information or access.

`docs/report.tex` (primary РПЗ) and `docs/REPORT.md` cite these IDs next to
material claims. A claim absent from this table is not a verified result.

| ID | Утверждение | Status | Evidence / scope |
|---|---|---|---|
| C-001 | Новый CPU solver содержит HLLD, MUSCL, staggered CT и SSP-RK2/Хойна | implemented | `src/kernels/*`, `src/MhdAmr.cpp` |
| C-002 | JSON inputs задают problem/grid/AMR/BC/time/scheme/output | implemented | `src/Config.cpp`, `inputs/*.json` |
| C-003 | Standalone Brio–Wu run: finite, 488 шагов, `max\|divB\|=0` | measured | `mhd2d_verify briowu`; регрессия `standalone.briowu` |
| C-004 | Standalone Alfvén `N=32`: `max\|divB\|=3.7e-13` | measured | `mhd2d_verify alfven 32`; регрессия `standalone.alfven32` |
| C-005 | AMR глобально консервативна | not established | нет gas reflux; фаза T07 |
| C-006 | Дискретный `divB` гарантирован для всей AMR-иерархии | not established | равномерный CT проверен; AMR — T07 |
| C-007 | Существует исполняемый GPU-путь | not established | горячие циклы — `amrex::LoopOnCpu`; фаза T12 |
| C-008 | H1: AMReX быстрее/точнее legacy | hypothesis → **partially supported** | см. C-020; равные resolution/cost/error кампании — T08+ |
| C-009 | Исторические результаты ВКР воспроизведены | not established | T03 = FAILED (`legacy_vkr` даёт NaN); показана согласованность порядка/диапазона, не побитовое воспроизведение — см. C-018/C-019 |
| C-010 | Экстремумы Brio–Wu классифицированы | not established | ветвь выбрана (compound/non-regular, D-003); независимый эталон + разметка фронтов — T06 |
| C-011 | Воспроизводимый CTest / run-manifest контур | verified | `ctest --preset release` = 9/9; `docs/REPRODUCIBILITY.md` |
| C-012 | `legacy_corrected` доводит все канонические задачи (кроме ОТ) до конца: finite, positive, conservative, div-free | measured | `docs/LEGACY_CORRECTED_T05_ENDTOEND.md`; `benchmarks/summary/legacy_corrected_endtoend.json`; 4 карточки, `quality_gate = pass` |
| C-013 | Исторический `legacy_vkr` расходится (отрицательные ρ,p) на Brio–Wu / CP-Alfvén / петле поля | verified | `docs/LEGACY_VKR_T03.md`; `scripts/run_legacy_vkr.py` quality gate |
| C-014 | `hlld_flux` имеет положительность-сохраняющий откат на HLL со счётчиком; на канонических тестах = 0 | implemented + verified | `src/kernels/Hlld.H`; `kernel.numerics::test_positivity_fallback`; `hlld_fallbacks` в диагностике всех прогонов |
| C-015 | Набор модульных/регрессионных тестов ядра AMReX | verified | `tests/kernel_numerics.cpp` (`kernel.numerics`); conversions, fast_speed, HLLD consistency/branches/degeneracies, 20k finiteness sweep, limiter, corner EMF |
| C-016 | Наблюдаемый порядок циркулярной альфвеновской волны ≥ 1.8 (N=64→128: 1.94 по относительной L1 B⊥; 1.92 по L1(sum4) — норма, которую проверяет гейт) | measured | `standalone.alfven_order`; `docs/figures/data/alfven_conv.dat` |
| C-017 | SSP-RK2 — второй порядок (полином усиления `1+a+a²/2`, порядок ОДУ 2.004) | verified | `kernel.numerics::test_ssprk2_order` |
| C-018 | AMReX ОТ/rotor совпадают с Авдеевой–Лукиным / ВКР по диапазонам (~2%); срез ОТ `y=0.3125` совпадает по признакам | measured (qualitative) | `docs/LITERATURE_VALIDATION.md`; `benchmarks/summary/rotor_slice_comparison.json`; `docs/figures/data/{ot_slice,rotor_diag_*}.dat` |
| C-019 | `legacy_corrected` петля поля: рост `divB` совпадает с ВКР §4.1.4 (до множителя ~2) | measured | `docs/LITERATURE_VALIDATION.md`; 1.4e-17→4.8e-16 vs 5.9e-17→2.0e-15 |
| C-020 | Brio–Wu 1D, равный шаг: AMReX N3 ~15× меньше `L1` и 113–222× быстрее legacy | provisional | `docs/CROSS_SOLVER_BRIOWU_1D.md`; `benchmarks/summary/briowu_1d_comparison.json`; эталон — сошедшийся AMReX N3 N=2048, не независимый (T06) |
| C-021 | `legacy_corrected` не доводит вихрь Орзага–Танга до конца при N=128 | verified | страж положительности при `t≈0.4314`; та же точка при CFL 0.5 и 0.2; `docs/LITERATURE_VALIDATION.md` |
| C-022 | На истинно 1D-сетке `HLLD_flux_corrected` (legacy) и `hlld_flux` (AMReX N0) совпадают ~6–10% на Brio–Wu и Dai–Woodward; обе строго монотонны | measured | `scripts/{run_legacy_1d_riemann.py,legacy_1d_riemann.cpp,compare_1d_riemann.py}`; `benchmarks/summary/{briowu_1d_flux,dai_woodward_1d}_comparison.json`; `docs/CROSS_SOLVER_COMPARISON.md §1` |
| C-023 | Петля поля: `legacy_corrected` сохраняет ≈18% E_B за два прохода, `mhd2d_amrex` ≈84% (N=128) → 0.91 (N=256) | measured | `mhd2d_verify loop`; `benchmarks/raw/legacy_corrected/magnetic_loop_scaled_128x64_cmp/`; `docs/figures/data/loop_eb_amrex.dat`; `docs/CROSS_SOLVER_COMPARISON.md §5` |
| C-024 | Тест Dai–Woodward добавлен в `mhd2d_verify` (`dw1d`) и legacy 1D-драйвер; обе схемы дают допустимое монотонное решение, диапазоны совпадают с ВКР рис. 7 | implemented + verified | `standalone.dai_woodward` (CTest); `docs/CROSS_SOLVER_COMPARISON.md §1` |
