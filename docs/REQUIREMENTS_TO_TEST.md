# Requirements-to-test register

Стартовая матрица T00, обновлено 2026-09-03. Статус `planned` означает, что
критерий принят, но ещё не доказан автоматическим тестом.

| ID | Требование | Проверка / будущий артефакт | Статус | Фаза |
|---|---|---|---|---|
| ENG-001 | clean configure/build | пресеты `release`, `cpu-release`, `cpu-debug`, `mpi-release`, `profile`, `hdf5-release` конфигурируются; `cuda-release` объявлен, но на этой машине не собирается (нет GPU) | **implemented (T02+T14)** | T02 |
| ENG-002 | reproducible regression | CTest 17/17: ядра, standalone, конфиг, канонические задачи, `amr.conservation`, `arch.kernel_purity`, manifest repeat. Проверено и под `cpu-debug` (AMReX assertions + bound check) | **implemented (T02+T04+T07)** | T02 |
| CFG-001 | strict JSON schema | valid and intentionally invalid JSON CTests | implemented (T02) | T02 |
| NUM-001 | primitive/conservative states | `kernel.numerics::test_conversions` (2000-sample round trip + `pressure_from_cons`) | **implemented (T04)** | T04 |
| NUM-002 | HLLD consistency and degeneracies | `kernel.numerics`: `F(q,q)=F_phys` (3000), supersonic branches, `Bn=0`, Brio--Wu, RD, 20000-sample finiteness sweep | **implemented (T04)** | T04 |
| NUM-003 | positivity observable and fallback | `hll_flux` fallback + validity predicate in `Hlld.H`; три раздельных счётчика `hlld_fallbacks` / `floor_events` (NEW-003, внутри `cons_to_prim`) / `nonpositive_cells`; `kernel.numerics::test_positivity_fallback`; гладкие задачи требуют нулей, МГД-взрыв даёт 16/48/0 | **implemented (T04)** | T04 |
| NUM-004 | uniform CT divergence | `divb: max_abs=… normalized=…` печатается в конце каждого прогона; порог `<=1e-12` проверяют `canonical.*` и `amr.conservation` | **implemented (T04)** | T04 |
| NUM-005 | SSP-RK2 coefficients and convergence | `kernel.numerics::test_ssprk2_order` (amplification `1+a+a²/2`, ODE order 2.004); `standalone.alfven_order` (observed 1.92 ≥ 1.8) | **implemented (T04)** | T04 |
| AMR-001 | gas conservation across coarse/fine | `YAFluxRegister` на каждую пару уровней; CTest `amr.conservation` (падает и при отсутствии стыка уровней, и если при `reflux=false` дефект не проявляется). Статическая иерархия: дрейф ρ 6.9e-5 → 6.3e-15, E 1.5e-4 → 7.9e-16 | **implemented (T07)** | T07 |
| AMR-002 | CT coarse/fine consistency | нормированный `div B` на иерархии проверяется в `amr.conservation` и `canonical.constant_state` (двухуровневая сетка) | **implemented (T07)** | T07 |
| LEG-001 | immutable `legacy_vkr` provenance | аннотированный тег `legacy_vkr/9d0f60e` + SHA-256 архива `afb51261…`, сверяются раннером на каждом прогоне | **implemented**; сами прогоны `legacy_vkr` по-прежнему дают NaN (T03 FAILED, ожидаемо) | T03 |
| LEG-002 | `legacy_corrected` physically admissible end-to-end | 5 карт доведены до конца (добавлен вращающийся цилиндр 128², манифест `benchmarks/raw/legacy_corrected/rotor_128/`); вихрь Орзага--Танга при N=128 не завершается, причина установлена (катастрофическое сокращение при β→0), см. `CROSS_SOLVER_COMPARISON.md` | **implemented (T05 CPU slice)** | T05 |
| BW-001 | Brio--Wu extrema diagnosis | независимый эталон (KT, без решателя Римана), версионированная разметка фронтов, метрики L1/L2/Linf + TV excess + over/undershoot + положение и ширина фронта | **implemented (T06)**: `docs/BRIOWU_T06.md`, гейт `briowu.independent_reference` | T06 |
| CMP-001 | legacy vs AMReX quality comparison | `scripts/compare_briowu_1d.py`, `benchmarks/summary/briowu_1d_comparison.json`, `docs/CROSS_SOLVER_BRIOWU_1D.md` (equal-spacing + equal-cost, provisional reference) | **partial (T08 Level A, Brio--Wu only)** | T08 |
| PERF-001 | benchmark provenance and repeats | JSON manifests, raw timings, median/spread | diagnostic timings only; pinned campaign absent | T09 |
| GPU-001 | CUDA parity | CUDA CTest suite and profiles | absent | T12 |
| RPT-001 | report claim evidence | [CLAIM_TO_EVIDENCE.md](CLAIM_TO_EVIDENCE.md) | planned | T01 |
