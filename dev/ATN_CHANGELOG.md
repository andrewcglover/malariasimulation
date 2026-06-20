# ATN Fork — Model Change Ledger

**Purpose:** crash-resilient, authoritative record of every model-affecting change made in this fork
since it diverged from upstream `mrc-ide/malariasimulation`. Primarily useful as a systematic audit
tool while diagnosing the relative ATN-vs-CFP impact discrepancy between this fork and the original
`malariasimple` deterministic model.

**Source of truth:** `git log` and `git show` on `atn-dev` since merge-base `5f63ccd`
(39 commits total; 11 model-affecting, listed below; remainder are plots, tests, site-file download,
session notes, RcppExports regen, gitignore, and pure documentation scaffolding).

**Status flags:** ✓ correct/intended · ⚠ suspect / open question · 🐛 was a bug, fixed in later commit

---

## Quick-reference table

| # | Commit | Date | Change (summary) | ATN impact direction |
|---|--------|------|-----------------|----------------------|
| 1 | `f989920` | 2026-06-15 | Add ATN params to `get_parameters()`, no-ATN defaults | Neutral (off by default) |
| 2 | `affd62b` | 2026-06-15 | Replace `Sm/Pm/Im` deque-EIP with 2-D Erlang `Sv/Ev/Iv`; `compute_atn_kernels()`; EIR=`sum(Iv)` | **Structural** |
| 3 | `e4ff574` | 2026-06-15 | Fix Erlang equilibrium init (ODE steady-state formula) | Correctness |
| 4 | `027ee7b` | 2026-06-16 | Derive `Lambda00 = foim·(1−B_max_post)`; remove `Lambda00sf`; `use_bompard` default→TRUE | Slightly **stronger** |
| 5a | `7e87ba0` | 2026-06-16 | **kappa `1/deltaq` → `1.0`** — much faster compartment drain | Materially **weaker** vs malsimple |
| 5b | `7e87ba0` | 2026-06-16 | Midpoint `i−0.5` → `i−1.5` (s=0.5 for first exposed compartment) | Slightly **stronger** |
| 6 | `69d9f1c` | 2026-06-16 | Add `atn_window`; derive `kappa = deltaq/atn_window`; `s = (i−1.5)·(atn_window/deltaq)` | Neutral at defaults; generalises #5 |
| 7 | `ce4f684` | 2026-06-17 | `phi_bednets[[species]]` indexing fix (was length-N vector crash) | Correctness |
| 8 | `0fd3c47` | 2026-06-18 | Site-sourced net retention; resistance-projected future `dn0/rn/gamman`; drop `pyr_atn dn0_atn` override | Context / Pyr-ATN trajectory |
| 9 | `08efe54` | 2026-06-18 | `lambda_atn` auto-derive = `1/retention` (was hardcoded 0) | **Weaker over time** (coverage wanes) |
| 10 | `d9cec67` | 2026-06-19 | Introduce `contact_factor`; `av_da = a·delta_atn·(sn+rnm)/sn` | Stronger (but `/sn` was buggy 🐛) |
| 11 | `48add43` | 2026-06-19 | Bound `contact_factor = (sn+rnm)/(1−rnm)` | **Stronger** (~+32%); fixes ordering |

---

## Per-change detail

### #1 — Add ATN parameters (`f989920`, 2026-06-15)

**Files:** `R/parameters.R`

**What:** Added all ATN parameters to `get_parameters()` with no-ATN defaults so existing call
signatures stay valid. Key defaults: `deltaq=1L`, `spor_len=10L`, `n_atn=1L`, `p_atn=0`,
`B_max_post=0`, `Lambda00sf=1` (removed later, see #4), `use_bompard=FALSE` (changed in #4).

**Impact:** None when ATN is off (`p_atn=0`). Sets the parameter surface the rest of the port builds on.

**Status:** ✓

---

### #2 — 2-D Sv/Ev/Iv Erlang ODE (`affd62b`, 2026-06-15)

**Files:** `src/adult_mosquito_eqs.h`, `src/adult_mosquito_eqs.cpp`, `R/biting_process.R`,
`R/compartmental.R`, `R/mosquito_biology.R`, `R/render.R`, `R/variables.R`

**What:** The structural core of the ATN port. Replaces:
- the legacy 3-compartment `Sm/Pm/Im` model with a `std::deque`-based discrete EIP delay
- with a 2-D `Sv[deltaqp1] / Ev[deltaqp1, spor_len] / Iv[deltaqp1]` Erlang ODE system

Key sub-changes:
- **C++ state vector** grows from 3 adult slots to `deltaqp1*(2+spor_len)` adult slots; new
  index helpers `sv_idx / ev_idx / iv_idx`; `kappa` hardcoded as `1.0/deltaq` (bug — see #5).
- **`compute_atn_kernels()`** added to `R/biting_process.R` — computes `delta_atn`, `Lambda_i`,
  `rho_i`, `B_post`, `dn_atn` per timestep and passes them to C++.
- **EIR read-out** changed from `solver_states[Im_idx]` to `sum(Iv block)`.
- **`ADULT_ODE_INDICES`** replaced by `sv/ev/iv_block_indices()` and `make_adult_ode_indices()`.
- Initial midpoint: `s = i - 0.5` (1-indexed, matching v3 at the time; bug — see #5).
- Initial kappa: `1.0/deltaq` (bug — see #5).

**Impact:** Structural. The ATN mechanism is now live; baseline (`delta_atn=0`) approximates stock
malariasimulation to tolerance (Erlang EIP converges to single-delay as `spor_len→∞`).

**Status:** ✓ structural correctness; sub-bugs in kappa + midpoint fixed in #5.

---

### #3 — Erlang equilibrium init (`e4ff574`, 2026-06-15)

**Files:** `R/mosquito_biology.R`

**What:** Replaced the legacy equilibrium init (which used the deque's `exp(-mu*dem)` incubation
survival) with the correct Erlang ODE steady-state:
```
E1_eq  = foim * n_Sv / (rho + mu)
Ev[j]  = E1_eq * (rho/(rho+mu))^(j-1)        j = 1..spor_len
n_Iv   = rho * E1_eq * (rho/(rho+mu))^(spor_len-1) / mu
```
Exposed rows are zeroed; ATN switches on at `t0_atn`. Also fixed `equilibrium_total_M()` to use
Erlang survival `(rho/(rho+mu))^spor_len` rather than the single-delay `exp(-mu*dem)`.

**Impact:** Init correctness — prevents an initial transient at t=0 in ATN-off runs.

**Status:** ✓

---

### #4 — `Lambda00` derivation; `use_bompard` default (`027ee7b`, 2026-06-16)

**Files:** `R/biting_process.R`, `R/parameters.R`

**What:**
```diff
- Lambda00 <- Lambda * parameters$Lambda00sf   # free scale-factor parameter
+ Lambda00 <- Lambda * (1 - parameters$B_max_post)   # tied to the fitted b_max
```
Removed the free `Lambda00sf` parameter. `Lambda00` — the pre-infection FOI floor at maximum
blocking — is now derived from the same `B_max_post` used for post-infection blocking, ensuring
the two drug effects share a single fitted parameter and cannot drift out of sync across posterior
draws.

Also changed `use_bompard` default from `FALSE` → `TRUE` (Bompard TRA→field-TBA transform active
by default; consistent with how Mali runs are configured).

**Impact:** Slightly **stronger** ATN — the new `Lambda00 = foim*(1−B_max_post)` with
`B_max_post ≈ 0.99999` gives a lower floor (less residual FOI at full blocking) than
`Lambda00sf = 1` (= no blocking at all, the old default). In practice this only matters when ATN
drug params are set; the no-ATN default `B_max_post=0` gives `Lambda00=foim` (identical to no
blocking), keeping the ATN-off path unchanged.

**Status:** ✓

---

### #5 — Kappa and exposure midpoint fixes (`7e87ba0`, 2026-06-16)

**Files:** `R/biting_process.R`, `R/parameters.R`, `src/adult_mosquito_eqs.h`,
`src/adult_mosquito_eqs.cpp`

This commit bundled **two opposing bug fixes**. Understanding them separately is essential.

#### 5a — kappa: `1.0/deltaq` → `1.0` (day⁻¹)

**Before:** `kappa` was hardcoded in the C++ constructor as `1.0 / deltaq`. With `deltaq = 10`,
`kappa = 0.1 day⁻¹` → mean residence **~10 days per exposed compartment**, **~100 days total**
in the drugged chain.

**After:** `kappa = 1.0` (now an explicit parameter passed from R, defaulting to 1.0) → ~1
day/compartment, ~10 days total — consistent with the `atn_window = 10` day exposure window.

**Note — this is NOT the same as v3/malsimple.** The reference (`…v3.R:484`) still uses:
```r
kappa <- 1 / deltaq    # = 0.1 with deltaq=10
```
v3/malsimple therefore has an **internal inconsistency**: its Hill time-axis spans s = 0.5→10.5
days (compartment midpoints), but its *physical* mean compartment residence is ~10 days, giving
~100 days total drugged persistence. The port corrected this to be consistent.

**Impact direction — the leading hypothesis for weaker ATN:**
- v3/malsimple: drug effect persists ~100 days in the exposed chain.
- Port: drug effect persists ~10 days.
- A mosquito in the drugged chain for ~100 days is **~10× more likely to survive to become
  infectious while still drug-affected** than one that drains in ~10 days. This is likely the
  **dominant reason ATN looks weaker in the port** relative to malsimple.

**Open question ⚠:** Was malsimple's 100-day persistence the *intended* biological effect (drug
lasts ~3 months in exposed mosquitoes) or was it itself an unintentional consequence of
`kappa = 1/deltaq`? The answer determines whether the port or malsimple is "correct" here.
The Hill time-axis and the compartment structure were designed to span ~10 days; the long
persistence appears inconsistent with that design, suggesting it was unintentional. But this
needs a deliberate decision before concluding the port's ~10-day persistence is right.

#### 5b — Exposure midpoint: `s = i − 0.5` → `s = i − 1.5`

**Before:** `s <- seq_len(deltaqp1) - 0.5` → first exposed compartment (i=2) got `s = 1.5`.
This matched v3's `(i - 0.5)` formula literally but was a **one-day over-ageing** of the
just-dosed cohort: a mosquito entering compartment 2 has been exposed for ≈ 0–1 days, so its
midpoint should be 0.5, not 1.5.

**After:** `s[1]=0; s[i] = (i-1.5)` → first exposed compartment at `s = 0.5`. Correct midpoint.

**Impact:** The Hill kernel is *decreasing* in `s` (higher `s` = weaker drug). So the fix
**increases** blocking for every exposed compartment — slightly **stronger** ATN. Opposes #5a.

**Note:** v3's `(i − 0.5)` therefore contains the over-ageing too; the port corrected both
models' shared latent offset. **Status: ✓ (port is correct; v3 still has the offset).**

---

### #6 — Add `atn_window`; generalise `kappa` derivation (`69d9f1c`, 2026-06-16)

**Files:** `R/biting_process.R`, `R/parameters.R`, `R/compartmental.R`

**What:** Introduced `atn_window` (default 10 days) as the total time-since-exposure window
spanned by the exposed compartments. Derived:
```r
kappa <- deltaq / atn_window       # at compartmental.R:77
s     <- (seq_len(deltaqp1) - 1.5) * (atn_window / deltaq)   # midpoint in days
```
The explicit `kappa` parameter from #5 was absorbed into this formula (removed from `parameters`).

**Impact at defaults (`atn_window = 10, deltaq = 10`):**
- `kappa = 10/10 = 1.0` — identical to #5; no change.
- `s spacing = 10/10 = 1 day` — identical to #5; no change.

**Why:** Decouples the resolution (number of compartments, `deltaq`) from the physical window
width (`atn_window`). Allows running higher `deltaq` for numerical convergence without
accidentally widening the drug-effect window. Also enables convergence testing.

**Status:** ✓

---

### #7 — `phi_bednets` species indexing fix (`ce4f684`, 2026-06-17)

**Files:** `R/biting_process.R`

**What:** `compute_atn_kernels()` previously used `parameters$phi_bednets` bare — a length-N
per-species vector after `set_species()`. Even `0 * phi_bednets` gives a length-N vector; Rcpp's
`as<double>()` rejects it with *"Expecting a single value: [extent=N]"*, crashing all Mali
parallel workers. Fixed by adding a `species` integer argument; all call sites updated.

**Impact:** Correctness — `delta_atn` is now correctly per-species (`phi_bednets[[species]]` is a
scalar). In single-species test runs this was a hard crash. In multi-species (Mali) it crashed all
36 parallel workers.

**Status:** ✓ (was a crash bug, not a modelling error in ATN strength)

---

### #8 — Site-sourced retention; resistance-projected future nets (`0fd3c47`, 2026-06-18)

**Files:** `dev/mali_projection_run.R` only (pipeline, not the package itself)

**What:**
- `retention_time` now read from `unique(site_obj$interventions$mean_retention)` (~2014 d for MLI)
  rather than a hardcoded 588 d.
- `build_future_schedule(arm, region)` looks up projected pyrethroid resistance from
  `site_obj$vectors$pyrethroid_resistance` (2000–2050 per-region table) for each distribution
  year and calls `med_net(pars, res_year)`. So `dn0/rn/gamman` now vary across the future window
  in step with rising resistance (e.g. Mopti 2025≈0.82 → 2031≈0.92).
- Dropped the `pyr_atn dn0_atn` override. Pyrethroid mortality for Pyr-ATN flows through the
  `dn0/rn` in the net schedule (resistance-projected); `dn0_atn` represents only the
  antimalarial's extra mortality (default 0) — avoids double-counting.

**Impact:** Mali pipeline context. Longer retention (2014 d vs 588 d) means ATN coverage decays
more slowly between campaigns. Rising resistance weakens pyrethroid efficacy in future years.
No effect on non-Mali single-species test runs.

**Status:** ✓

---

### #9 — `lambda_atn` auto-derive from net retention (`08efe54`, 2026-06-18)

**Files:** `R/parameters.R`, `R/vector_control_parameters.R`, `R/biting_process.R`

**What:** Changed `lambda_atn` default from `0` (no coverage waning) to `NULL` (auto-derive
sentinel). `set_bednets()` now sets `lambda_atn = 1/bednet_retention` automatically when
`lambda_atn` is NULL. In `compute_atn_kernels()`, a NULL guard falls back to `lambda_atn = 0`
when no `set_bednets` call has been made.

Rationale: `set_bednets` models net loss as stochastic `Exponential(mean = retention)` via
`log_uniform`. `exp(-lambda_atn * t)` with `lambda_atn = 1/retention` is the deterministic
mean-field analogue — same exponential decay, same retention parameter.

In Mali: `retention_time ≈ 2014 d` → `lambda_atn ≈ 0.000497 day⁻¹`. ATN/Pyr-ATN coverage now
decays between campaigns. Previously (`lambda_atn = 0`) coverage was flat after each distribution.

**Impact:** **Weaker ATN over time** (between campaigns, coverage decays). Before this change ATN
coverage was implicitly permanent once distributed — an overestimate between campaign rounds.

**Status:** ✓ (correct modelling of net loss; the change makes ATN *less* overoptimistic)

---

### #10 — `contact_factor` introduced; initial formula (`d9cec67`, 2026-06-19)

**Files:** `R/biting_process.R`

**What:** Previously `av_da = a * delta_atn`. A mosquito physically repelled by the net
(probability `rnm`, the untreated-net barrier floor) still touches the net and picks up the
drug; only chemical excito-repellency (`rn − rnm`, the insecticide-driven component) prevents
contact entirely. The fix:
```r
av_da = a * delta_atn * contact_factor
contact_factor = (sn + rnm) / sn     # initial formula — BUGGY (see #11)
```
Sourced per ATN event from the bednet schedule (matching `t0_atn → bednet_timesteps`); falls back
to 1 when no `set_bednets` call.

**Impact:** **Stronger** — barrier-repelled mosquitoes now flow into `Sv_exposed`. Only `av_da`
changes; `delta_atn` fraction (for `Sv→Ev` infection split), `foim`, EIR, and total mosquito
density are unaffected.

**Bug in initial formula 🐛:** `sn = 1 − rn − dn → 0` for a strongly insecticidal net; `/sn`
made `contact_factor ≈ 2.2` for fresh Pyr-ATN vs ≈1.32 for non-insecticidal ATN — **backwards**
(Pyr-ATN showed more exposure than ATN despite pyrethroid repellency). Fixed in #11.

---

### #11 — Bounded `contact_factor` denominator (`48add43`, 2026-06-19)

**Files:** `R/biting_process.R`

**What:** Fixed the `/sn` blow-up from #10:
```r
# Before (buggy):
cf_each <- (sn_e + rnm_e) / sn_e

# After (bounded):
cf_each <- pmax(sn_e + rnm_e, 0) / pmax(1 - rnm_e, 1e-6)
```
Denominator `(1 − rnm)` is the untreated-net floor (~0.76) — fixed, bounded, never collapses.
Numerator `(sn + rnm) = 1 − rn_chem − dn` excludes pyrethroid-killed mosquitoes (cannot transmit).

**Behaviour:**
- Non-insecticidal ATN (`rn0 = rnm`, `dn0 = 0`): `1/(1 − rnm) ≈ 1.316` (constant, same as buggy
  version for this arm — main ATN result unchanged).
- Pyr-ATN fresh: `(1 − rn_chem − dn)/(1 − rnm) < 1/(1 − rnm)` — correctly **less** exposure than
  bare ATN; decays toward `1/(1 − rnm)` as insecticide wanes.
- No-ATN arms (`cfp`, `none`): `delta_atn = 0 ⇒ av_da = 0` regardless — bit-identical.

**Impact:** **Stronger** than feed-only (correctly includes barrier-repelled), but **weaker** than
the buggy #10 formula for Pyr-ATN. Restores Pyr-ATN < ATN exposure ordering.

**Status:** ✓

---

## Net-impact synthesis

Tallying the direction of each change relative to what malsimple/v3 would compute:

| Factor | Direction | Magnitude |
|--------|-----------|-----------|
| #5a kappa 10× faster drain | **Weaker** vs malsimple | **Large** — leading hypothesis |
| #9 `lambda_atn` coverage waning | **Weaker** over time | Moderate (between campaigns) |
| #5b midpoint anchor corrected | Stronger | Small (1-day shift on Hill) |
| #4 `Lambda00` from `B_max_post` | Stronger | Small |
| #11 `contact_factor` ~1.32 | Stronger | Moderate (~32% more exposure) |

The **kappa** change (#5a) is the dominant candidate. Both remaining weakeners (#5a and #9)
represent genuine corrections to over-optimistic assumptions in v3/malsimple:

- **kappa:** malsimple's ~100-day drugged residence may have been an unintentional consequence of
  `kappa = 1/deltaq` given the Hill axis spans only ~10 days.
- **lambda_atn:** malsimple had no net-loss decay on ATN coverage; permanent coverage is unrealistic.

The three strengtheners in the port at least partially offset, so the remaining net weakness is
primarily attributable to the kappa difference.

---

## Open questions

1. **kappa / drug persistence (⚠ most important):** In malsimple, `kappa = 1/deltaq = 0.1` with
   `deltaq = 10` gives ~100 days total drugged residence; in the port, `kappa = 1.0` gives ~10 days.
   The Hill kernel spans s = 0 to `atn_window = 10` days in both — a 10× mismatch in v3. Is the
   intended biological persistence ~10 days or ~100 days? Options:
   - **Keep port as-is** (10 days): consistent, matches the Hill axis design.
   - **Revert to `kappa = 1/deltaq`**: replicates malsimple; adds the same inconsistency but
     enables a direct like-for-like comparison.
   - **New `kappa` from biology**: set `atn_window` to the actual intended drug-persistence window
     (e.g. ~30 days) and derive `kappa = deltaq/atn_window` accordingly. This is the cleanest
     approach if there is a biologically-grounded estimate.

2. **`lambda_atn` in malsimple:** Did the original malsimple simulations use `lambda_atn = 0`
   (permanent coverage) or some decay? If zero, this accounts for additional ATN strength there.

3. **Contact factor — v3 analogue:** v3/malsimple has no `contact_factor` (exposure = feed rate
   only). The port's ~+32% from barrier-repelled mosquitoes is a genuine mechanistic addition with
   no counterpart in malsimple. Is this intended as a port improvement, or should it be disabled for
   a direct like-for-like comparison run?

---

## Baseline-correctness note

The ATN-off path (`delta_atn = 0`, set by `p_atn = 0`) leaves all ATN compartments inert:
- `Lambda_i[q] = Lambda` for all q (no blocking).
- `av_da = 0` (no exposure flux); only baseline compartment (q=0) is populated.
- EIR = `sum(Iv)` = `Iv[0]`, equal to the stock `Im` (to Erlang tolerance).

The Erlang EIP converges to the single-delay result as `spor_len → ∞`; with `spor_len = 10`
the match is very close but not bit-identical. This is expected and pre-existing behaviour.

---

*Generated 2026-06-20. Source: `git log/show` on `atn-dev` since `5f63ccd`. Update when new
model-affecting commits land on `atn-dev`.*
