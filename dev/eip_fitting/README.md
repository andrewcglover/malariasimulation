# EIP-slowdown kernel fitting (archive)

Stan models + R drivers used to fit the ATN **EIP-suppression** kernel, collected here for
provenance. These fit the sporozoite-positivity assay data (`SPZ_pooled_summary_APR26.csv`,
paired control vs 200 mg arms, dissection days 10/13) and produce the posterior consumed by
malariasimulation (`dev/atn_params/eip_fit_hill_result.rds`).

## Files and origin

| File | Kernel | Origin repo |
|---|---|---|
| `eip_fit_hill.stan` | **Hill** recovery `ρ(s)=ρ−(ρ−ρ0)·H^η/(H^η+s^η)` — **the version used in production** | `malariasimple_ATNs/dev` (active) |
| `fit_eip_hill.R` | driver for the Hill fit | `malariasimple_ATNs/dev` (active) |
| `eip_fit.stan` | **Exponential** recovery `ρ(s)=ρ−(ρ−ρ0)·exp(−ζs)` — earlier variant, superseded | `malariasimple_ATNs_190526backup/dev` (archived; removed from the active repo) |
| `fit_eip.R` | driver for the exponential fit | `malariasimple_ATNs_190526backup/dev` (archived) |
| `priors_atn_main.R` | shared prior definitions + `write_priors_json()`; `source()`d by both drivers | `malariasimple_ATNs/dev` (active) |

## Notes

- Both variants share everything except the recovery kernel: the exponential uses one parameter
  `ζ`; the Hill uses two (`s_half`, `nH`). Same Erlang(R=10, 1) completion, same paired-binomial
  likelihood with a per-row (observation-level) logit offset, same 7-day baseline-EIP floor.
- The **Hill** variant was chosen for production (`use_eip_hill = TRUE`); the exponential is kept
  here for reference/comparison.
- The assay data CSV is **not** included here (lives in `malariasimple_ATNs/dev/exp_data/`); to
  re-run a fit, set `data_path` accordingly and `source("priors_atn_main.R")` from this folder.
