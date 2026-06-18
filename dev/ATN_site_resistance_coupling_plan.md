# Plan: site-sourced future net params + full feeding-cycle ATN repellency coupling

## Context

Auditing `dev/mali_projection_run.R` and the ATN mosquito kernel surfaced three gaps. Confirmed findings:

1. **Retention & future net efficacy are not site-sourced.** `retention_time` is hardcoded at
   588 days (line 33), ignoring the site file. Future net params (`dn0`/`rn`/`gamman`) come from a
   resistance-grid lookup keyed to a single *historical-median* resistance per region (line 169),
   not the site file's projected resistance. Genuine projections exist in
   `site_obj$vectors$pyrethroid_resistance` (per-region, per-year, **2000–2050**, rising trend) —
   the earlier check had looked at `interventions` (which only carries 2024 forward). Past nets
   already read per-row site `dn0`/`rn0`/`gamman` correctly; only the future half and retention need fixing.

2. **ATN repellency ignores pyrethroid resistance.** `delta_atn = p_atn * phi_bednets[[species]] * Q_t`
   (`R/biting_process.R:312`) has no repellency/survival term. Pyrethroid resistance reaches the ATN
   mechanism only indirectly via the global biting rate `a`. The reference's intended fix is the
   commented 8-category `w/z` feeding model (`dev/reference/malariasimple_aitn_deterministic_v3.R:925–994`):
   the ITN+ATN category carries pyrethroid-driven `w = 1 - phi_bednets + phi_bednets*s_itn` and
   `z = phi_bednets*r_itn`. Chosen approach: **full feeding-cycle coupling**.

3. **Species-specific contact probability** — already satisfied: `compute_atn_kernels` indexes
   `phi_bednets[[species]]`, and Part B's coupling adds per-species `rn`/`dn` (`bednet_rn[matches, species]`),
   reinforcing it. No standalone work; verify only.

Outcome: Mali future projections driven by site-file retention and projected resistance, and an ATN
model where (especially Pyr-ATN) net repellency/mortality is consistently resistance- and
species-dependent through the same feeding-cycle accounting as the ITN model — while remaining
bit-compatible with upstream when ATN is off (`p_atn = 0`).

## Step 0 — Commit & push current state first (prerequisite)

Snapshot the working tree before any new edits. ATN work goes to **`atn-dev`**, not `master`.

- Currently uncommitted: `R/compartmental.R` (prior-session rendering rewrite),
  `dev/mali_projection_run.R`, `dev/mali_single_region_test.R`, and untracked
  `dev/reference/InterventionExpansion.R`.
- `git add -A` → commit (message describing the rendering aggregation + Mali plot/test changes) →
  `git push -u origin atn-dev` (creates/updates the remote `atn-dev` branch; the local branch
  currently tracks `origin/master`, so push explicitly to `atn-dev`, do **not** push to `master`).
- Confirm the push succeeded before starting Part A.

## Part A — Site-source retention + project future resistance (`dev/mali_projection_run.R`)

**A1. Retention from site file, with manual override.**
- Add near the top settings block: `retention_override <- NULL  # numeric days to override site value`.
- After `site_obj` is loaded, read the site value: `site_retention <- unique(site_obj$interventions$mean_retention)`
  (the `site` package mandates a single value; `site:::add_itns` does the same `unique(...)`). Guard
  `length(site_retention) == 1`.
- `retention_time <- if (!is.null(retention_override)) retention_override else site_retention`.
- Leave the CD coverage math (lines 43–45) and the `set_bednets(retention = retention_time)` call
  unchanged — they now consume the site value. **Flag in a comment**: site `mean_retention` (~2014 d)
  is far longer than the old 588 d and materially changes CD top-up coverage; this is the intended
  site-driven behaviour, override available if needed.

**A2. Future net efficacy from projected resistance.**
- Build a per-region year→resistance lookup from `site_obj$vectors$pyrethroid_resistance`
  (columns `name_1`, `year`, `pyrethroid_resistance`; covers 2000–2050).
- Change `build_future_schedule(arm, res)` → `build_future_schedule(arm, region)`. For each net
  distribution timestep on the future grid, map `year = start_year + floor(timestep / 365)`, look up
  that region/year resistance, and call `med_net(pars, res_year)` (median posterior draw conditional
  on that year's resistance). `dn0`/`rn`/`gamman` become per-distribution vectors aligned to the grid.
- Update the call site in `build_params` (line 184) to pass `region`. The single-median `res`
  (line 169) is no longer used for the future schedule.
- `med_net()` is unchanged (already takes a scalar resistance and returns medians); it is now called
  per distribution year instead of once.

**A3. `dn0_atn` (Pyr-ATN ATN-kernel mortality).**
- With Part B routing pyrethroid kill/repellency through the ITN-side kernel (`sn`/`rn`, now
  resistance-projected via A2), the ATN kernel's `dn0_atn` must represent **only the antimalarial's**
  extra mortality, not the pyrethroid's, to avoid double-counting. Drop the pyrethroid-derived
  `atn_overrides$dn0_atn <- med_net(only_pars, res)$dn0` assignment (line 196) for `pyr_atn`; leave
  `dn0_atn` at its antimalarial default (0 unless a true antimalarial-mortality estimate is supplied).

## Part B — Full feeding-cycle coupling for ATN repellency (`R/biting_process.R`)

Goal: make ATN drug-contact and Pyr-ATN pyrethroid effects flow through the same resistance-/species-
dependent feeding-cycle accounting (`w`/`z`) the ITN model already uses, per the reference's 8-category
expansion (`malariasimple_aitn_deterministic_v3.R:925–994`).

**B1. Thread the net kernel into `compute_atn_kernels`.**
- At the biting-loop call site (`R/biting_process.R:204`), these are already in scope per species/timestep:
  `f`, `a`, `W` (`average_p_successful`), `Z` (`average_p_repelled`), `Q0`, and `p_bitten` (from
  `prob_bitten`, carrying per-individual `prob_bitten_survives`/`prob_repelled`, i.e. the
  resistance-decayed `sn = 1 - rn - dn` and `rn` from `prob_survives_bednets`/`prob_repelled_bednets`).
- Extend `compute_atn_kernels(timestep, parameters, foim, species)` to also accept the feeding-cycle
  summaries (`W`, `Z`, and/or the net-contact terms) so `delta_atn` and the ATN exposure are derived
  from the same `w`/`z` weights as the ITN biting rate, rather than the raw `phi_bednets * Q_t`
  independent of the pyrethroid kernel.

**B2. Reformulate `delta_atn` per the reference 8-category model.**
- Implement the ITN+ATN feeding category: contact/exposure weighted by pyrethroid survival
  (`w = 1 - phi_bednets + phi_bednets*s_itn`) and repellency (`z = phi_bednets*r_itn`), with ATN
  coverage `Q_atn_t` coupled to ITN coverage on shared nets (Pyr-ATN) rather than independent. The
  scalar `delta_atn` passed to the ODE encodes the population-average antimalarial-contact-per-feeding
  consistent with `a`/`f` — no per-individual structure is needed in the deterministic ODE.
- Keep the change **R-side** in `compute_atn_kernels`, feeding results through the existing
  `adult_mosquito_model_update` arguments (`delta_atn`, `dn_atn`, `av_da = a * delta_atn`). Only touch
  C++ (`src/adult_mosquito_eqs.*`) if a needed combined term genuinely cannot be folded into the
  existing arguments — prefer to avoid C++ changes.
- The exact `w`/`z` algebra (whether `delta_atn`'s formula changes vs. only the `a`/anthropophagy
  coupling) will be confirmed against the reference lines during implementation; the hard constraint
  is the baseline-invariance check below.

**B3. Update call site and tests.**
- Pass the new feeding-cycle args at `R/biting_process.R:204` (production uses `s_i`).
- Update test call sites that invoke `compute_atn_kernels` directly with the new signature:
  `tests/testthat/test-atn-mosquito.R:33,92`, `tests/testthat/test-compartmental.R:52` (pass species `1L`
  and representative `W`/`Z`).

## Concern 3 — species-specific contact (verify only)

No code change. After Part B, confirm `phi_bednets[[species]]` and `bednet_rn[matches, species]` /
`bednet_dn0[matches, species]` make the ATN kernel fully per-species. The Mali `set_bednets` call
(`mali_projection_run.R:242–244`) currently replicates `dn0`/`rn` across species via
`matrix(rep(...), ncol = n_sp)`; note this in passing — true per-species future efficacy would require
per-species resistance, out of scope unless the site file provides it.

## Verification

1. **Baseline invariance (must pass):** with `p_atn = 0`, the fork reproduces upstream EIR/prevalence
   to tolerance. Run `tests/testthat/` — `test-atn-mosquito.R`, `test-compartmental.R`, and the
   stock-comparison baseline test must stay green after the signature change.
2. **Rebuild after any C++/signature change:** `Rcpp::compileAttributes()` then `devtools::load_all()`.
3. **Pyr-ATN resistance sensitivity (new behaviour):** in `dev/mali_single_region_test.R` (Mopti,
   `mali_test_mode`), confirm the `pyr_atn` arm now responds to projected resistance — higher resistance
   reduces pyrethroid killing/repellency and shifts ATN drug-contact accordingly — vs the `atn` and
   `cfp` arms. Inspect the aggregated `Sv/Ev/Iv_exposed_*` rendering columns.
4. **Site-sourcing checks:** assert `retention_time` equals the site `mean_retention` (or the override),
   and that future-grid `dn0`/`rn`/`gamman` vary across the 2025–2030 window in step with
   `site_obj$vectors$pyrethroid_resistance` for the region.
5. Run the single-region 4-arm harness end-to-end before any full 36-job parallel sweep.

## Files

- `dev/mali_projection_run.R` — retention override+site source (A1), future-resistance projection (A2),
  drop pyrethroid `dn0_atn` override (A3).
- `R/biting_process.R` — `compute_atn_kernels` signature + `delta_atn` reformulation (B1, B2); call
  site at line 204 (B3).
- `tests/testthat/test-atn-mosquito.R`, `tests/testthat/test-compartmental.R` — updated call signature (B3).
- `src/adult_mosquito_eqs.cpp` / `.h` — only if a combined term cannot be folded into existing
  `adult_mosquito_model_update` args (avoid if possible).
- Reference (read-only): `dev/reference/malariasimple_aitn_deterministic_v3.R:925–994` (8-category w/z).
