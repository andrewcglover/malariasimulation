# Session notes — 2026-06-18

## Where we got to

Full Mali admin-1 sweep (`dev/mali_projection_run.R`) completed successfully: **32/36 jobs saved**
to `dev/outputs/mali_projection_results.rds`. Results are ready for `mali_projection_plots.R`.

---

## Outstanding issue: Bamako ODE solver failure

**All 4 arms (none / cfp / atn / pyr_atn) fail for Bamako.**

Error (confirmed by the new per-solver label added this session):

> Solver error: too much work in the **adult mosquito** ODE solver.

### What this tells us

- All 4 arms fail — including `none` (no nets at all). This is **not ATN-related**. The adult
  ODE for Bamako cannot complete a single day's integration within `ode_max_steps = 1e7` sub-steps
  regardless of intervention arm.
- The corrected error label confirms it is the **adult solver** (not aquatic, despite the old
  hardcoded message that said "aquatic life stage model" — that label was wrong and has been fixed).
- Root cause is almost certainly Bamako's **seasonality profile**: a very sharp rainfall peak drives
  an extreme carrying-capacity transient that makes the adult ODE stiff on that day.
- The adult solver state is 27 states (non-ATN, `deltaq=1, spor_len=10`) vs 3 in upstream
  malariasimulation. This fork's Erlang chain carries more states through the same solver budget.

### Immediate next step

Two options to try **in order**:

1. **Raise `ode_max_steps` to `1e8` in `build_params()`** (`dev/mali_projection_run.R:240`).
   Quick — just change `1e7` → `1e8` and rerun. If Bamako completes, done.
2. **Investigate Bamako's seasonality parameters** — check whether Bamako has an unusually peaked
   rainfall profile in the site file (`site_obj$seasonality`) that creates a pathological carrying-
   capacity spike on a specific day. Possible fix: small smoothing or floor on carrying capacity.
   Only needed if 1e8 is still insufficient.

---

## Changes made this session (all committed and pushed to `origin/atn-dev`)

### Commits (newest first)

| Hash | Description |
|------|-------------|
| `08efe54` | ATN net-retention decay + correct ODE solver-failure attribution |
| `034e42f` | Mali sweep: cap workers at 8, raise ode_max_steps, capture per-job failures |

### `08efe54` — detail

**A. ODE solver error label fixed (C++)**
- `src/solver.h` / `src/solver.cpp`: shared `Observer` now accepts a `model_name` string.
  Error message says "adult mosquito" or "aquatic mosquito larval" correctly.
- `src/adult_mosquito_eqs.cpp` / `src/aquatic_mosquito_eqs.cpp`: pass the label at construction.
- Error dump now also prints `n_states` — adult >> aquatic (27/135 vs 3) as a cross-check.

**C. `lambda_atn` auto-derive from net retention (R, package-level)**
- `R/parameters.R`: `lambda_atn` default changed from `0` to `NULL` (auto-derive sentinel).
- `R/vector_control_parameters.R` (`set_bednets`): when `lambda_atn` is NULL, sets
  `lambda_atn = 1 / bednet_retention` automatically. For logistic retention: warns and maps
  `1 / bednet_logistic_half_life`.
- `R/biting_process.R` (`compute_atn_kernels`): NULL guard → falls back to 0 if no `set_bednets`.
- Effect in Mali pipeline: `build_params()` calls `set_bednets(retention = retention_time)`
  → `lambda_atn = 1/2013.6` ≈ 0.000497 days⁻¹. ATN/Pyr-ATN coverage now wanes between
  campaigns (previously pinned flat). No change to `mali_projection_run.R` needed.
- `tests/testthat/test-atn-mosquito.R` §12e: 4 new tests; all 30 ATN tests pass.

### `034e42f` — detail

- `N_CORES` capped at 8 (`min(8L, detectCores()-1)`).
- `ode_max_steps = 1e7` override in `build_params()` (band-aid; Bamako still fails — see above).
- `run_one` wrapped in `tryCatch`: returns `list(.__error__=TRUE, region, arm, message)` on error.
- After sweep: failures printed + written to `dev/outputs/mali_projection_results_failures.csv`.
- `dev/outputs/` added to `.gitignore`.

---

## The three ATN decay components (confirmed correct in code)

| # | Component | Mechanism | Status |
|---|-----------|-----------|--------|
| 1 | Repellency / insecticide decay | `gamman` → `prob_bitten` per-individual via `net_time` | ✅ Always active |
| 2 | Antimalarial drug-activity decay | `gamma_atn` → per-event `Lambda0`/`rho0`/`dn_atn`, coverage-weighted | ✅ Active (set in Mali pipeline) |
| 3 | Net retention / loss of nets | `lambda_atn` → `exp(-lambda_atn*t)` on `Q_t` per event | ✅ **Now active** (auto-derived from `bednet_retention`) |

**Key insight on component 3:** `set_bednets` models net loss as stochastic
`Exponential(mean = retention)` via `log_uniform` (`R/utils.R:28`). The ATN coverage decay
`exp(-lambda_atn * t)` with `lambda_atn = 1/retention` is the deterministic mean-field analogue —
same exponential, differing only in stochasticity.

---

## Parallel execution setup (confirmed working)

- 8 workers via `parallel::makeCluster` (PSOCK, Windows).
- Workers load pre-compiled `.dll` via `pkgload::load_all(compile=FALSE)` — so **must**
  `devtools::load_all()` (compile) in RStudio before launching sweep after any C++ change.
- Log is silent during the sweep (`parLapplyLB` blocks until all jobs finish; worker stdout/stderr
  not forwarded). Only the main process messages appear at start and end.
- Output: `dev/outputs/mali_projection_results.rds` + `_failures.csv` if any failures.
- Safe to run single-region test in RStudio console during sweep (different output files,
  DLL read-only load is not locked).

---

## To pick up tomorrow

1. Fix Bamako: try `ode_max_steps = 1e8` in `build_params()`, rerun sweep.
2. If Bamako resolves: run `mali_projection_plots.R` on the full results.
3. If Bamako still fails: investigate `site_obj$seasonality` for Bamako.
