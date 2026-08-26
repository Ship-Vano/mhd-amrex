# Requirements-to-test register

Стартовая матрица T00. Статус `planned` означает, что критерий принят, но
ещё не доказан автоматическим тестом.

| ID | Требование | Проверка / будущий артефакт | Статус на T00 | Фаза |
|---|---|---|---|---|
| ENG-001 | clean configure/build | CMake preset `release` | implemented in T02 | T02 |
| ENG-002 | reproducible regression | CTest: kernel unit, standalone Brio--Wu/Alfvén, config validation and manifest repeat | implemented in T02; numerical scope remains limited | T02 |
| CFG-001 | strict JSON schema | valid and intentionally invalid JSON CTests | implemented in T02 | T02 |
| NUM-001 | primitive/conservative states | unit tests including invalid states | planned | T04 |
| NUM-002 | HLLD consistency and degeneracies | unit/regression tests with branch coverage | planned | T04 |
| NUM-003 | positivity observable and fallback | counters + failure diagnostics | planned | T04 |
| NUM-004 | uniform CT divergence | `L_inf` and normalized norm in CTest artifact | manual diagnostic only | T04 |
| NUM-005 | SSP-RK2 coefficients and convergence | Butcher and smooth-wave tests | planned | T04 |
| AMR-001 | gas conservation across coarse/fine | periodic conservation test with flux register | known unmet | T07 |
| AMR-002 | CT coarse/fine consistency | AMR `divB` regression | manual diagnostic only | T07 |
| LEG-001 | immutable `legacy_vkr` provenance | tag, checksum, three canonical runs | blocked by D-001/D-002 | T03 |
| BW-001 | Brio--Wu extrema diagnosis | fixed reference + metrics + annotations | blocked by D-003 | T06 |
| PERF-001 | benchmark provenance and repeats | JSON manifests, raw timings, median/spread | absent | T09 |
| GPU-001 | CUDA parity | CUDA CTest suite and profiles | absent | T12 |
| RPT-001 | report claim evidence | [CLAIM_TO_EVIDENCE.md](CLAIM_TO_EVIDENCE.md) | planned | T01 |
