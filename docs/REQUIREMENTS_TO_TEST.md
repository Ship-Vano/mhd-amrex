# Requirements-to-test register

Стартовая матрица T00, обновлено 2026-08-31. Статус `planned` означает, что
критерий принят, но ещё не доказан автоматическим тестом.

| ID | Требование | Проверка / будущий артефакт | Статус | Фаза |
|---|---|---|---|---|
| ENG-001 | clean configure/build | CMake preset `release` | implemented (T02) | T02 |
| ENG-002 | reproducible regression | CTest 9/9: kernel unit + numerics, standalone Brio--Wu/Alfvén/alfven_order, config validation, manifest repeat | implemented (T02+T04) | T02 |
| CFG-001 | strict JSON schema | valid and intentionally invalid JSON CTests | implemented (T02) | T02 |
| NUM-001 | primitive/conservative states | `kernel.numerics::test_conversions` (2000-sample round trip + `pressure_from_cons`) | **implemented (T04)** | T04 |
| NUM-002 | HLLD consistency and degeneracies | `kernel.numerics`: `F(q,q)=F_phys` (3000), supersonic branches, `Bn=0`, Brio--Wu, RD, 20000-sample finiteness sweep | **implemented (T04)** | T04 |
| NUM-003 | positivity observable and fallback | `hll_flux` fallback + validity predicate in `Hlld.H`; `hlld_fallbacks` / `nonpositive_cells` counters in both drivers; `kernel.numerics::test_positivity_fallback`; smooth cases assert 0 | **implemented (T04)** | T04 |
| NUM-004 | uniform CT divergence | `L_inf` and normalized norm in a CTest artifact | printed only; artifact still planned | T04 |
| NUM-005 | SSP-RK2 coefficients and convergence | `kernel.numerics::test_ssprk2_order` (amplification `1+a+a²/2`, ODE order 2.004); `standalone.alfven_order` (observed 1.92 ≥ 1.8) | **implemented (T04)** | T04 |
| AMR-001 | gas conservation across coarse/fine | periodic conservation test with flux register | known unmet — no gas reflux | T07 |
| AMR-002 | CT coarse/fine consistency | AMR `divB` regression | manual diagnostic only | T07 |
| LEG-001 | immutable `legacy_vkr` provenance | tag, checksum, three canonical runs | `legacy_vkr` runs are provenance-clean but reproduce NaN (T03 FAILED, expected) | T03 |
| LEG-002 | `legacy_corrected` physically admissible end-to-end | 4 canonical cards run to completion, finite/positive/conservative/div-free; `docs/LEGACY_CORRECTED_T05_ENDTOEND.md` | **implemented (T05 CPU slice)** | T05 |
| BW-001 | Brio--Wu extrema diagnosis | fixed reference + metrics + annotations | branch chosen (D-003 = compound); provisional 1-D comparison done; independent reference still pending | T06 |
| CMP-001 | legacy vs AMReX quality comparison | `scripts/compare_briowu_1d.py`, `benchmarks/summary/briowu_1d_comparison.json`, `docs/CROSS_SOLVER_BRIOWU_1D.md` (equal-spacing + equal-cost, provisional reference) | **partial (T08 Level A, Brio--Wu only)** | T08 |
| PERF-001 | benchmark provenance and repeats | JSON manifests, raw timings, median/spread | diagnostic timings only; pinned campaign absent | T09 |
| GPU-001 | CUDA parity | CUDA CTest suite and profiles | absent | T12 |
| RPT-001 | report claim evidence | [CLAIM_TO_EVIDENCE.md](CLAIM_TO_EVIDENCE.md) | planned | T01 |
