# Plan — ATN drug-exposure should reflect *net contact*, not *successful feeds*

## Context

While inspecting single-site Mali projection runs we noticed total mosquito counts were
**higher under the future-ATN arm than under the no-future-nets ("none") arm**, for all species.
We traced that anomaly to pre-existing leMenach/Griffin feeding-cycle logic on the `set_bednets`
side (a repellent-only net, `dn0=0`/`rn=0.24`, lowers the per-day death rate `mu` via the feeding
rate `f`). **That anomaly is explicitly out of scope here** (see §"Out of scope").

This plan addresses a *separate, genuine* modelling refinement we surfaced in that discussion:
the ATN **drug-exposure** rate currently equals `a * delta_atn`, where `a` is the *successful-feed*
rate. But a mosquito picks up the antimalarial whenever it **contacts the net**, not only when it
feeds. A mosquito that is physically repelled by the net (the historical untreated-net floor,
`rnm`) still *touched the net → got dosed*. Only **chemical excito-repellency** (`rn − rnm`, the
insecticide-induced part that decays toward 0) keeps a mosquito off the net entirely.

So today's model **under-doses**: it excludes barrier-repelled mosquitoes from the exposed pool.
The fix makes ATNs more realistically effective at blocking transmission, with **zero change** to
mosquito density, human infection, EIR, or the existing baseline.

## Design decision (settled with user)

**Approach 2 ("true contact rate"), Refinement 2 ("scale the rate only").**

- The repellency relevant to *drug exposure* is the **chemical** part only, `rn − rnm`. The physical
  floor `rnm` represents net contact, so barrier-repelled mosquitoes **are** dosed.
  - Non-insecticidal ATN (`rn0 = rnm`, `dn=0`): chemical repellency ≡ 0 → maximal exposure.
  - Pyr-ATN (`rn0 > rnm`, `dn>0`): chemical repellency starts at `rn0 − rnm` and decays to 0 with
    net age — exposure recovers as the insecticide wanes. Falls out of the existing decay
    `rn(dt) = (rn0 − rnm)·exp(−dt/gamman) + rnm` automatically; `rnm` is the dynamic floor (no
    hardcoded 0.24, no new "0.24" constant).

- **Change only the exposure *rate*** (`a → contact rate` inside `av_da`). **Keep the `delta_atn`
  *fraction* and the FOI terms (`Lambda_i`, `Lambda0_t`) tied to feeding** — because *dosing needs
  contact* but *infection needs a blood meal*. This is the key correctness point: the extra
  (barrier-repelled, non-feeding) mosquitoes flow into `Sv_exposed` only, **never** into
  `Ev_exposed` (you cannot be infected without feeding). Verified against the ODE structure:
  `Sv→Ev[1]` rate `= delta_atn·Lambda0_t·Svtot` stays unchanged; `Sv→Sv[1]` absorbs the entire
  increase; the split still balances (it is not held constant — total exposure outflow rises).

### The contact factor

Per net-user encounter: `P(feed&survive) = sn = 1 − rn − dn`, `P(barrier-repelled) = rnm`,
`P(killed) = dn`. Contacts that survive to transmit = feeds + barrier-repelled, so

```
contact_factor = (sn + rnm) / sn = (1 − rn − dn + rnm) / (1 − rn − dn)
               = (1 − rn_chem − dn) / (1 − rn − dn),   rn_chem = rn − rnm
```

- Non-insecticidal ATN: `(1 − 0 − 0)/(1 − 0.24 − 0) = 1/0.76 ≈ 1.316` (constant).
- Pyr-ATN: varies with net age via `rn(dt)`, `dn(dt)`; → `1/(1−rnm)` as the net ages.
- Killed mosquitoes (`dn`) are excluded (already removed via elevated `mu`/reduced `a` on the
  `set_bednets` side — do **not** re-count them).
- **Numerical guard:** floor `sn` (e.g. `max(sn, 1e-6)`) to avoid divide-by-zero if `rn+dn → 1`.

Exposure rate becomes `av_da = a · delta_atn · contact_factor`. With `delta_atn = 0` (ATN-off),
`av_da = 0` regardless — **baseline untouched**.

## Implementation

All changes are **R-side**; the C++ ODE (`src/adult_mosquito_eqs.cpp`) is **unchanged** — it still
receives `av_da` and the bare `delta_atn` as two separate arguments; we only change what we pass
for `av_da`.

1. **`compute_atn_kernels()` — `R/biting_process.R:278`** (primary change)
   - Compute a coverage-weighted `contact_factor` over the ATN distribution events, reusing the
     existing per-event age/decay machinery (`Q_each`, `age`, `bednet_decay`-style terms) that
     already produces `Lambda0_t`/`rho0_t`/`dn_atn`.
   - **Source the net's `rn`/`rnm`/`dn0`/`gamman` for each ATN event** by matching `t0_atn` against
     the bednet schedule (`parameters$bednet_timesteps`) and indexing
     `parameters$bednet_rn[ , species]` / `bednet_rnm` / `bednet_dn0` / `bednet_gamman`
     (per-species matrices; use `species`, mirroring the existing `phi_bednets[[species]]` rule —
     CLAUDE.md §10). `t0_atn` aligns with the future portion of the schedule by construction
     (`mali_projection_run.R:214-215, 248-258`).
   - **Fallback:** if no `set_bednets` call has been made (net arrays NULL — same situation the
     `lambda_atn` NULL guard handles), set `contact_factor = 1` (no adjustment). Keeps direct
     unit-test calls and ATN-off paths valid.
   - Add `contact_factor` (scalar) to the returned list.

2. **Per-step feed-in — `R/biting_process.R:209`**
   - Change the `av_da` argument from `a * kernels$delta_atn` to
     `a * kernels$delta_atn * kernels$contact_factor`.
   - **Leave `R/biting_process.R:210` (`kernels$delta_atn`) exactly as is** — the fraction `da`
     stays feed-based for the `(1−da)·Lambda_i` / `da·Lambda0_t` infection-splitting terms.

3. **No new `get_parameters()` keys** if sourcing from the bednet schedule (preferred). If matching
   proves brittle, the fallback is explicit per-event ATN net vectors (`rn_atn`/`rnm_atn`/…)
   passed alongside `Q0_atn`/`t0_atn`; document the choice in CLAUDE.md §10 either way.

## What does NOT change (guard rails)

- **C++** (`adult_mosquito_eqs.cpp`/`.h`), state layout, RcppExports — untouched.
- **`a` used for `foim`** (`biting_process.R:175`, human→mosquito infection into baseline `Sv[0]`)
  and **EIR** (`calculate_eir`, `infectious * a`) — untouched (recomputed from the original `a`).
- **`mu`, `f`, emergence, aquatic/larval model, the entire `set_bednets`/biting side** — untouched
  (the mosquito-count anomaly therefore persists; that is intended here).
- **Baseline / ATN-off** runs — bit-for-bit unchanged (`delta_atn = 0 ⇒ av_da = 0`).

## Verification

1. **Unit test (new), `tests/testthat/test-atn-mosquito.R` §12d-style** — call `compute_atn_kernels()`
   directly and assert:
   - non-insecticidal ATN (`rn0 = rnm`, `dn0 = 0`) → `contact_factor == 1/(1−rnm)`;
   - higher `dn`/`rn_chem` (fresh Pyr-ATN) → smaller surviving-contact term, → `1/(1−rnm)` as age ↑;
   - no `set_bednets` → `contact_factor == 1` (fallback);
   - `delta_atn` (the fraction) is **unchanged** by this work.
2. **Baseline regression** — existing baseline test still passes (ATN-off identical to upstream).
3. **Single-region sim — `dev/mali_single_region_test.R` (Mopti)** — run `atn` and `pyr_atn`:
   - confirm `Sv_exposed`/`Ev_exposed`/`Iv_exposed` counts **rise** vs current (more dosing);
   - confirm `Ev_exposed` inflow is *not* disproportionately inflated (the extra goes to
     `Sv_exposed`), and PfPR/EIR drop modestly more under ATN (stronger blocking);
   - confirm `none`-arm and total mosquito **density** is unchanged vs current (sanity: the anomaly
     is untouched, not "fixed").
   Run via `devtools::load_all()` then the test script (`options(mali_test_mode = TRUE)`).

## Out of scope (flag for future)

- The **mosquito-count anomaly** (ATN > none) is the leMenach/Griffin repellency-lowers-`mu`
  behaviour on the `set_bednets` side. The user has decided to leave the full leMenach/Griffin
  death-rate logic as-is for now and revisit it in a future piece of work.
- A fuller "contact rate" that also recomputes the feeding cycle (Approach 1) was rejected — it
  perturbs `f`/gonotrophic timing and risks other dynamics.

## Follow-ups after approval

- Copy this plan to `dev/` (e.g. `dev/atn_contact_exposure_plan.md`) for project-tracked reference.
- Update **CLAUDE.md §3** ("Repellency / pyrethroid-resistance coupling") and §10 to record this as
  a deliberate *revision* of the earlier "all repellency flows through `a`" decision: barrier
  repellency (`rnm`) is net contact ⇒ dosed; only chemical repellency (`rn − rnm`) suppresses ATN
  exposure. Note `contact_factor` and that it scales the `av_da` rate only.
- Commit to `atn-dev` (per memory: ATN work lands on `atn-dev`, merge to `master` at milestones).
