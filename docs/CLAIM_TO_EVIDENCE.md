# Claim-to-evidence matrix for the NIRS draft

Версия: T01, 2026-08-26. Статусы имеют строгий смысл:

- `implemented` — код присутствует; это не численное доказательство;
- `diagnostic` — единичное наблюдение, не benchmark и не полная verification;
- `planned` — критерий и будущая проверка определены;
- `blocked` — нужна внешняя информация или доступ.

| ID | Утверждение в РПЗ | Status | Evidence / scope |
|---|---|---|---|
| C-001 | Новый CPU solver содержит HLLD, MUSCL, staggered CT и SSP-RK2/Heun | implemented | source tree at `40ebbae`; [BASELINE_MANIFEST.md](BASELINE_MANIFEST.md) |
| C-002 | JSON inputs задают problem/grid/AMR/BC/time/scheme/output | implemented | `src/Config.cpp`, `inputs/*.json` at `40ebbae` |
| C-003 | Standalone Brio--Wu smoke run is finite, 488 steps, printed `max|divB|=0` | diagnostic | command and output in [BASELINE_MANIFEST.md](BASELINE_MANIFEST.md); one local run |
| C-004 | Standalone Alfvén `N=32` printed `max|divB|=3.706e-13` | diagnostic | command and output in [BASELINE_MANIFEST.md](BASELINE_MANIFEST.md); one local run |
| C-005 | AMR is globally conservative | not established | gas reflux absent; T07 required |
| C-006 | Discrete divergence is guaranteed for the whole AMR hierarchy | not established | uniform diagnostic is insufficient; T04/T07 required |
| C-007 | GPU execution path exists | not established | hot path uses `amrex::LoopOnCpu`; T12 required |
| C-008 | AMReX is faster/better than legacy | hypothesis | H1 requires equal-resolution/cost/error campaigns, T08+ |
| C-009 | Legacy VCR results are reproduced | not established | external archive lacks exact run provenance; T03 blocked by D-001/D-002 |
| C-010 | Brio--Wu oscillations are classified | not established | reference branch not selected; T06 blocked by D-003 |
| C-011 | Reproducible CTest/run-manifest foundation | planned | T02 |

The two report sources, `REPORT.md` and `report.tex`, cite these IDs next to
all material claims. A claim absent from this table must not be interpreted as
a verified result.
