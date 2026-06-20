# CLAUDE.md — malariasimulation fork: add ATN (antimalarial-treated net) functionality

> Place this file in the **root of the malariasimulation fork** (`you/malariasimulation`, branch `atn-dev`).
> Claude Code reads it automatically each session.
>
> **Maintenance (Claude Code):** keep this file current as you work. Whenever you confirm a source fact,
> locate a file / function / constructor signature, settle a design decision, or hit a non-obvious gotcha,
> update the relevant section concisely and replace any "TO CONFIRM" note with the finding — so future
> sessions inherit it. Keep edits tight; this is a working reference, not prose.

## 0. One-line goal

Fork `mrc-ide/malariasimulation` and port the **antimalarial-treated net (ATN)** mechanism
from the author's `malariasimple` deterministic model into malariasimulation's **deterministic
adult-mosquito ODE**, keeping the result an installable, backward-compatible R package that runs
on the Imperial DIDE cluster via `hipercow`.

## 1. Ground-truth source — READ THIS FIRST

The authoritative ATN model is the working `malariasimple` implementation (odin2/dust2,
discrete-time). Read it in full before editing anything:

- `malariasimple_aitn_deterministic_v3.R`  ← the current canonical version (1149 lines)
  - Mosquito / ATN block: **lines ~351–625** (states, kernels, dSv/dEv/dIv, totals)
  - ATN distribution-events block: **lines ~773–864**
  - Larval block (shared, unchanged): lines ~628–685

Keep a copy of this file inside the fork (e.g. `dev/reference/`) so Claude Code can always read it.
When in doubt about an equation, the v3 file is correct, not this summary.

## 2. What changes vs what stays

malariasimulation = **individual-based humans + deterministic compartmental ODE mosquitoes**
(aquatic `E/L/P` + adult `Sm/Pm/Im`, integrated in C++ via the adult/aquatic solver classes).

**Only the adult mosquito infection block changes.** Everything else — aquatic stages, larval
seasonality, the entire human IBM, immunity, biting/EIR machinery — stays as-is. The two models
share the same lineage (Griffin 2010 / White 2011 / Imperial), so the larval and human sides are
already equivalent.

## 3. The ATN mechanism (what to build)

Replace the 3-compartment adult model (`Sm/Pm/Im`) with a **2-D ATN-exposure structure**:

- `Sv[deltaqp1]`            susceptible, by ATN-exposure compartment
- `Ev[deltaqp1, spor_len]`  latent/sporogony, ATN-exposure × Erlang stage
- `Iv[deltaqp1]`            infectious, by ATN-exposure compartment

Indices:
- `deltaqp1 = deltaq + 1`. Compartment **1 = unexposed (baseline)**; **2..deltaqp1 = days since
  ATN exposure**. A conveyor of rate `kappa = 1/deltaq` carries mosquitoes from exposed
  compartments back toward baseline (drug effect wanes over `deltaq` days).
- `spor_len` = Erlang sub-stages approximating the EIP; baseline progression `rho = spor_len/delayMos`.
  **This Erlang chain is the single biggest structural change** — malariasimulation's single `Pm`
  latent compartment cannot represent time-since-infection, which `rho_i` and `B_post` (below) need.

**Exposure coupling:** `delta_atn = p_atn * phi_atn * Q_atn_t`, the probability a mosquito meets
the antimalarial on an attempted bite. A mosquito biting an ATN host (rate `av*delta_atn`) moves
into exposure compartment 2. `Q_atn_t` = total ATN coverage across distribution events.

**Repellency / pyrethroid-resistance coupling — DESIGN DECISION (revised 2026-06-19).**
The exposure rate is `av_da = a * delta_atn * contact_factor` (see `R/biting_process.R`,
`compute_atn_kernels`). The three terms are:

- **`a`** (human blood-meal rate, `R/biting_process.R:118`) carries the full feeding-cycle
  adjustment from `W`/`Z`, embedding per-individual `sn`/`rn`/`dn` via `prob_survives_bednets` /
  `prob_repelled_bednets`. Pyrethroid resistance ↑ ⇒ `rn`/`dn` ↓ ⇒ `a` ↑ ⇒ more ATN exposure.
- **`delta_atn`** (the *fraction* `p_atn * phi_bednets * Q_t`) carries no repellency term — it is
  the net-user coverage exposure probability, not a rate. The FOI-splitting terms (`Lambda_i`,
  `Lambda0_t`) remain tied to `a` (feeding) because *infection requires a blood meal*.
- **`contact_factor`** (new, 2026-06-19; corrected formula 2026-06-19) corrects the exposure *rate*
  for mosquitoes that **physically touch the net but do not feed**: barrier-repelled mosquitoes
  (prob `rnm`, the untreated-net floor) touch the net → pick up the drug, so they belong in the
  exposed pool. Only **chemical excito-repellency** (`rn − rnm`, the insecticide-driven part) keeps
  a mosquito off the net entirely. Formula (bounded, excludes pyrethroid-killed):
  `contact_factor = (sn + rnm) / (1 − rnm)`, where `sn = 1 − rn − dn`, `rn_chem = rn − rnm`.
  - Denominator `(1 − rnm)` is the untreated-net floor (~0.76) — always bounded, never collapses.
    *Previous formula* `(sn + rnm)/sn` was buggy: `/sn → 0` for insecticidal nets inflated Pyr-ATN
    exposure above ATN, inverting the correct ordering. Fixed to `/(1 − rnm)`.
  - Numerator `(sn + rnm) = 1 − rn_chem − dn` excludes pyrethroid-killed mosquitoes (`dn`); a
    dead mosquito cannot transmit.
  - Scaling the realized `a` (not a no-net `a0`) preserves IRS + historical-net coupling for free:
    IRS and historical non-ATN bednets already suppress `a` via `W`/`Z`; only the ATN net's own
    repellency/mortality is re-applied via `contact_surv`. No double-counting.
  - Non-insecticidal ATN (`rn0 = rnm`, `dn0 = 0`): `(1−rnm)/(1−rnm) = 1/(1−rnm)` (constant,
    same as before — main ATN result unchanged by the formula correction).
  - Pyr-ATN fresh net: `(1 − rn0 − dn0 + rnm)/(1 − rnm)`; decays to `1/(1 − rnm)` as insecticide
    wanes (`rn → rnm`, `dn → 0`). Correctly `< 1/(1−rnm)` (Pyr-ATN < ATN, ordering restored).
    Derived dynamically per event from the bednet schedule (`parameters$bednet_rn`/`rnm`/`dn0`/
    `gamman`, matched by `t0_atn`).
  - ATN-off (`delta_atn = 0`): `av_da = 0` regardless — **baseline untouched**.
  - No `set_bednets` call: falls back to `contact_factor = 1`.
  - **Only `av_da` is affected.** `foim` (→ `Sv[0]`), EIR (`calculate_eir`), `mu`/`f`, the
    aquatic model, and total mosquito density are all unchanged. `Sv→Ev` infection rates stay
    feed-based; the extra contacts from `contact_factor` flow into `Sv_exposed` only.

**Out-of-scope note:** the leMenach/Griffin death-rate formula `p1 = p1_0·W/(1 − Z·p1_0)` means
a repellent-only net (ATN, `dn0=0`, `rn=0.24`) lowers `mu` relative to no net, so total mosquito
density is slightly *higher* under ATN than no-nets. This is pre-existing biting-model behaviour,
not an ATN-port artifact; revisiting the full leMenach/Griffin death-rate logic is deferred.

The full per-net-category `av_mosq[i] = av*w[i]/wh` (v3 commented 8-category block,
lines ~925–994) is a larger structural change — not pursued.

**Coverage is ATN-only by design.** `Q_atn_t` is built solely from the ATN distribution events
(`Q0_atn`/`t0_atn`) in `compute_atn_kernels`, so it correctly **excludes** historical pyrethroid
ITNs still in circulation after the first ATN campaign — those suppress biting (via per-individual
`net_time` in `a`) but deliver no drug. Do **not** couple `delta_atn`'s coverage to the `W`/`Z`
net-using population: that population includes the historical ITN users and would mis-attribute drug
exposure to them during the ITN→ATN transition.

**Four drug effects** (all decay across exposure compartments via Hill kernels, optionally
Bompard TRA→field-TBA transformed; see v3 lines ~430–467):

1. **Pre-infection blocking** `Lambda_i` — reduces the human→mosquito FOI for exposed mosquitoes.
2. **EIP suppression** `rho_i` — slows sporogony for exposed mosquitoes (longer EIP ⇒ fewer reach `Iv`).
3. **Extra mortality** `dn_atn` — added death for exposed mosquitoes (the `(1 - dn_atn)` survival factors).
4. **Post-infection blocking** `B_post[j]` — an already-infected `Ev` mosquito, re-exposed, clears its
   infection and returns to `Sv[2]`; keyed to time-since-infection `t_post[j] = (j-0.5)*delayMos/spor_len`.

**Drug potency decay over calendar time:** each ATN distribution event's effect (`Lambda0`, `rho0`,
`dn`) relaxes toward baseline at rate `gamma_atn` since its `t0_atn`. Multiple events: `n_atn`
events with vector `t0_atn`, `Q0_atn`; coverage-weighted averaging + random proportional
replacement for overlapping coverage (v3 lines ~773–864).

**Coupling back to the rest of the model:**
- human → mosquito FOI (`foim`) feeds `Sv[1]` (unexposed).
- mosquito → human EIR uses `Ivtot = sum(Iv[])` in place of `Im`.

## 4. Implementation plan (route A — edit in place, keep it a package)

1. **Locate the adult mosquito ODE** in `src/` (the adult-mosquito solver class + its derivative
   function) and the R objects that construct it. Inspect `R/RcppExports.R` to find every C++↔R
   entry point (constructor signature, state read-out).
2. **Rewrite the adult derivative + state vector** to the `Sv/Ev/Iv` system above. Adult state grows
   from 3 to `deltaqp1*(2 + spor_len)`. Transcribe dSv/dEv/dIv from v3 lines ~538–603, but:
   - **drop the `dt` factors** — those are discrete-Euler hacks; the continuous adaptive solver
     integrates them natively. The `if (x + dx < 0) 0` guard from odin2/malariasimple is
     replaced by the `nn` lambda non-negativity floor in `create_eqs` (see §8).
3. **Plumb the new parameters** through `get_parameters()` (R) → the C++ constructor (Rcpp).
   Give every ATN parameter a **no-ATN default** (`Q0_atn = 0`, `n_atn = 1`, etc.) so existing call
   signatures and the `site`/`cali`/`scene`/`postie` ecosystem keep working unchanged.
4. **Equilibrium / initialisation** (the main gotcha): the enlarged state breaks malariasimulation's
   closed-form mosquito equilibrium. Use the v3 approach — seed the **baseline compartment (index 1)**
   at the no-ATN steady state, zero all exposed compartments, let ATN switch on at `t0_atn`. Confirm
   the baseline collapses exactly to stock `Sm/Pm/Im`.
5. **Decide where Hill/Bompard kernels live.** Simplest: compute `delta_atn`, `Lambda_i`, `rho_i`,
   `dn_atn`, `B_post` R-side each step (like ITN net-decay) and pass vectors into C++; the C++ side
   then only does the ODE arithmetic.
6. Run `Rcpp::compileAttributes()` + `devtools::load_all()` / `install()` after C++ interface changes.

## 5. Verification (do not skip)

- **Baseline test:** with ATN off, the fork must reproduce **stock malariasimulation** EIR and
  prevalence to tolerance. The baseline compartment must behave identically to `Sm/Pm/Im`.
- **ATN-on test:** with ATN on, compare against `malariasimple_aitn_deterministic_v3.R` under matched
  parameters (same `deltaq`, `spor_len`, kernels, coverage, `t0_atn`).
- Add a regression test under `tests/testthat/` for both.

## 6. Environment, build & deployment

- **Local:** Windows + RStudio + **Rtools** (matching the R version). Compile via
  `devtools::load_all()` / `devtools::install()`. Claude Code runs in RStudio's Terminal pane.
- **Isolation:** develop in a dedicated RStudio Project; use `devtools::load_all()` (never install the
  fork over the user's stable malariasimulation used in other projects). `renv` for full isolation +
  reproducible cluster runs is optional but recommended.
- **Cluster:** Imperial **DIDE** cluster via **hipercow** (drivers `dide-windows` and `dide-linux`).
  Provision the fork from GitHub: `pkgdepends` with `you/malariasimulation@atn-dev` (or a tag), or the
  `script` method (`remotes::install_github(...)`). `conan2` compiles it on the cluster. C++ is plain
  portable ODE arithmetic, so Windows or Linux nodes both work. Re-provision after each C++ change you
  want on the cluster.

## 7. Git workflow (Claude Code can run all of this)

- Origin = the user's fork; add `upstream` = `mrc-ide/malariasimulation`.
- `main` = stable, deployable line; `atn-dev` = active work.
- Merge `atn-dev → main` only when it compiles **and** passes the baseline test.
- Sync upstream by merging `upstream/main → atn-dev` first (test), then `→ main`. Conflicts will be
  small (edits are localised to the mosquito files).
- **Tag `main` at milestones** (e.g. `v0.1-atn`) so hipercow can provision exact, reproducible versions.

## 8. Watch-outs

- EIP-as-delay → Erlang chain is the structural crux; get the `rho`/`rho_i` stage rates right.
- Keep backward compatibility: no-ATN defaults must make the model bit-identical to upstream.
- Continuous ODE ≠ discrete Euler: remove `dt` discretisation factors from the dxdt expressions.
  **Exception — non-negativity floor (added 2026-06-20):** a `nn` lambda in `create_eqs`
  (`src/adult_mosquito_eqs.cpp`) clamps all adult compartment reads to `max(x[i], 0)`. This is
  required because high `a_tol` (set in the Mali pipeline to prevent near-zero thrashing)
  permits large steps that overshoot to slightly negative, triggering a NaN cascade via
  wrong-sign outflow terms. The floor is robustness, not speed (may marginally slow near-zero
  crossings). See `dev/bamako_solver_diagnosis.md`.
- Don't rename the package — it would break the `malariasimulation::`-calling ecosystem.
- **Dimension bounds:** require `deltaqp1 ≥ 2` (i.e. `deltaq ≥ 1`, since `kappa = 1/deltaq`),
  `spor_len ≥ 1`, and `n_atn ≥ 1`. For an ATN-*off* run use `n_atn = 1` with `Q0_atn = 0` (keep
  `deltaqp1 ≥ 2`, `delta_atn = 0`) — never collapse the dimension with `deltaq = 0`.
- **Empty-loop caveat:** the reference's `3:deltaqp1` and `2:spor_len` slices must become C-style loops
  (`for (i = 3; i <= deltaqp1; ++i)`) that are *naturally empty* at the minimum sizes — don't write code
  assuming `deltaqp1 ≥ 3` or `spor_len ≥ 2`. odin2 does this automatically; the hand-port must too.
- **Structural constants:** `deltaqp1`, `spor_len`, `n_atn` set array sizes — fixed at model
  construction (rebuild to change), not runtime knobs. Keep `deltaqp1 == deltaq + 1`; better, compute
  `deltaqp1 = deltaq + 1` internally and expose only `deltaq`.

## 9. Confirmed source facts (verified by reading the files, 11 Jun 2026)

The adult mosquito model was inspected directly. These supersede the "locate the…" guidance in §4
step 1 — the files and structure below are known, not to-be-discovered.

**Files (both in `src/`):**
- `adult_mosquito_eqs.cpp` — the derivative + the Rcpp-exported entry points.
- `adult_mosquito_eqs.h` — the state enum, struct, constructor declaration.

**Current structure (what's being replaced):**
- State is one flat vector: aquatic `E/L/P` at indices 0–2 (from `aquatic_mosquito_eqs.h`), adult at
  3–5 via `enum class AdultState : size_t {S = 3, E = 4, I = 5}`. `get_idx` just casts the enum to
  `size_t`, so enum values *are* the vector positions.
- **EIP is a discrete delay, not a compartment:** a `std::deque<double> lagged_incubating` of length
  `tau` in the `AdultMosquitoModel` struct. `adult_mosquito_model_update` pushes `susceptible*foim`
  and pops the front; the derivative moves `lagged_incubating.front() * exp(-mu*tau)` from E→I. This
  is exactly why the Erlang `Ev` chain is a genuine structural change (the deque can't carry per-stage
  `rho_i` / `B_post`).
- Aquatic coupling already sums the adult states: `total_M = S + E + I` (cpp lines ~29–32).

**Exported `[[Rcpp::export]]` functions (signatures change → RcppExports regenerates):**
`create_adult_mosquito_model`, `adult_mosquito_model_update`, `adult_mosquito_model_save_state`,
`adult_mosquito_model_restore_state`, `create_adult_solver`.

**Target index layout (proposed, keep aquatic at 0–2):**
- `Sv[q]   → 3 + q`                                   (q = 0..deltaq)
- `Ev[q,j] → 3 + deltaqp1 + q*spor_len + j`           (j = 0..spor_len-1)
- `Iv[q]   → 3 + deltaqp1*(1 + spor_len) + q`
- total state size = `3 + deltaqp1*(2 + spor_len)`.

**C++ touch-list:**
- `.h`: replace the `AdultState` enum with parameterised index helpers (need `deltaqp1`,`spor_len`);
  in the struct, remove `lagged_incubating` and `const double tau`, add `deltaq`/`deltaqp1`/
  `spor_len`, `kappa`, `rho`, and per-step ATN terms (`delta_atn`, `dn_atn`, vectors `Lambda_i`,
  `rho_i`, `B_post`); update the constructor decl; add a state-size accessor.
- `.cpp`: rewrite `create_eqs` to the `Sv/Ev/Iv` system; replace `total_M = S+E+I` with sums; drop the
  `incubation_survival`/lagged logic; widen `create_adult_mosquito_model` and
  `adult_mosquito_model_update`; `save_state`/`restore_state` likely collapse (no deque to checkpoint);
  `create_adult_solver` body unchanged (only its `init` is longer).

**R-side touch-points (confirmed):**
- `get_parameters()` — `R/parameters.R:341` — add ATN params + no-ATN defaults.
- **mosquito equilibrium / init — `set_equilibrium()` at `R/parameters.R:264`**: calls
  `malariaEquilibrium::human_equilibrium()` for the *human* state only; extracts `eq$FOIM` →
  `init_foim`. All mosquito-compartment arithmetic is malariasimulation's own:
  `parameterise_mosquito_equilibrium()` → `initial_mosquito_counts()` (`R/mosquito_biology.R:10`) →
  `create_adult_solver()`. Only `initial_mosquito_counts()` needs to emit the enlarged `Sv/Ev/Iv`
  init vector (whole equilibrium into baseline index 1 via `Ev_ratio = rho/(rho+mu)` &
  `Ev_norm_factor`, exposed rows zero). Do **not** fork or modify `malariaEquilibrium`. Because
  `malariaEquilibrium` — hence `init_foim` / total density — is unchanged, baseline init just
  *redistributes the same total*; but it reproduces stock only **to tolerance** (single-delay EIP
  vs Erlang; converges as `spor_len` grows), not bit-identically.
- **`ADULT_ODE_INDICES`** — `R/compartmental.R:1–2` — the R-side mirror of the C++ `AdultState`
  enum: `c(Sm = 4, Pm = 5, Im = 6)`. Must expand in step with the C++ enum when the state vector
  widens to `Sv/Ev/Iv`.
- **Per-step update feed-in** — `R/biting_process.R:204` — `adult_mosquito_model_update(model, mu, foim, Sm, f)` called each timestep; `foim` routes into `Sv[1]` (baseline only) after the rewrite.
- **EIR read-out** — `R/biting_process.R:262` — currently reads `solver_states[[ADULT_ODE_INDICES['Im']]]` (single index); after rewrite becomes a sum over the entire `Iv` block (`Ivtot = sum(Iv[1..deltaqp1])`).

**Compatibility strategy:** disaggregate internally, **sum only at the boundary** — `sum(Iv)` for the
EIR read-out and `sum(Sv)/sum(Ev)/sum(Iv)` for any `Sm/Pm/Im` outputs — and route `foim` specifically
into `Sv[1]` (exposed rows use `Lambda_i`, not baseline `foim`).

**Reference → package name map** (v3 symbol → malariasimulation equivalent):

| v3 / malariasimple | malariasimulation | Notes |
|---|---|---|
| `delayMos` | `dem` (`parameters$dem`) | EIP duration (days) |
| `mu` | `mum` (`parameters$mum`) | adult mosquito death rate |
| `Lambda` / `foim` in v3 | `foim` / `init_foim` | human→mosquito FOI; `init_foim = eq$FOIM` at equilibrium |
| `mv0` | `m` (local in `initial_mosquito_counts`) | total adult mosquito density |

## 10. Naming, collisions & parameter conventions

- **Only one shared namespace matters:** the `get_parameters()` named list. C++ struct members and
  locals (`Sv`, `Ev`, `Iv`, `kappa`, `rho`, `deltaq`) are scoped and cannot clash with anything else.
  R lists *silently tolerate duplicate names* (`list(Q0 = 1, Q0 = 2)` is valid), so a key collision
  won't error — it just misbehaves. Hence the check below.
- **Before adding ATN params:** enumerate the existing `get_parameters()` keys, then classify each
  reference quantity as (a) a NEW parameter to add, (b) an EXISTING one to reuse (`mu`, `tau`/`delayMos`,
  `phi_bednets`), or (c) a derived quantity. Check every new name against the existing list.
- **Convention:** suffix ATN-specific *parameter-list keys* with `_atn` (already started: `Q0_atn`,
  `t0_atn`, `p_atn`, `gamma_atn`); keep *C++ internals* short and unsuffixed (they're scoped, and the
  dense math reads better). Verbose names live in the param list where they document intent; tight names
  live in the C++.
- `Q0` already exists in malariasimulation (anthropophagy / human blood index); ATN coverage is
  `Q0_atn`, already distinct — no clash.
- **ATN distribution events:** single and multiple rounds use the *same* code path — `n_atn` with
  vectors `t0_atn`, `Q0_atn` (`n_atn = 1` for a single round). Proportional-replacement mixing assumes
  `t0_atn` is chronological and distinct.
- **`use_bompard` default is TRUE** (changed 2026-06-16); `use_eip_hill` default is also TRUE.
  Both only affect ATN-on runs; baseline tests are unaffected.
- **`Lambda00sf` was removed** (2026-06-16). `Lambda00` is now derived inside `compute_atn_kernels()`
  as `foim * (1 - B_max_post)` — pre- and post-infection blocking share the same `b_max`. There is no
  `Lambda00sf` parameter; passing it as an override will error.
- **`p_atn`** (default 0) must be set to a positive value (e.g. 0.9) when distributing ATNs;
  it multiplies `phi_bednets * Q_t` to give `delta_atn`. Omitting it silently disables the ATN mechanism.
- **`set_bednets` requires `rnm < rn` (strict)**. For non-insecticidal ATNs (`dn0 = 0`, `rn = 0.24`),
  set `rnm = rn - 1e-9`; `rnm = rn` is physically correct but fails the API check.
- **`phi_bednets` is per-species** (set by `set_species` to a length-N vector). `compute_atn_kernels`
  takes a `species` integer argument and must use `parameters$phi_bednets[[species]]` — NOT the bare
  `parameters$phi_bednets` vector. Even `0 * phi_bednets` yields a length-N vector, causing
  `Expecting a single value: [extent=N]` from Rcpp when passed as the scalar `delta_atn` argument.
  All call sites must pass the species index (production: `s_i`; tests: `1L`).
- **`lambda_atn` is now NULL by default (auto-derive sentinel, changed 2026-06-18).** `set_bednets`
  automatically sets `lambda_atn = 1/bednet_retention` from whatever retention value it received
  (manual or site-sourced) — so ATN drug coverage wanes at the same rate as net loss on the
  pyrethroid side. Rationale: `set_bednets` net loss is stochastic Exponential(mean=retention)
  via `log_uniform` (`R/utils.R:28`); the ATN decay `exp(-lambda_atn*t)` is its deterministic
  mean-field analogue. An explicit `lambda_atn` value (e.g. `0`) in `get_parameters(overrides=...)`
  is respected and NOT overwritten. For logistic retention, a warning is issued and
  `lambda_atn = 1/bednet_logistic_half_life` is used. `compute_atn_kernels` falls back to
  `lambda_atn=0` (no waning) if no `set_bednets` call has been made. Tests in §12e.
- **`contact_factor` sourcing (added 2026-06-19; formula corrected 2026-06-19).** `compute_atn_kernels`
  reads `rn`/`rnm`/`dn0`/`gamman` for each ATN event via `match(t0_atn, parameters$bednet_timesteps)`.
  In the Mali pipeline every ATN event matches exactly one bednet row (same `timesteps` vector,
  `mali_projection_run.R:248-258`). Formula: `(sn_e + rnm_e) / (1 − rnm_e)`, where
  `sn_e = 1 − rn_e − dn_e`, then coverage-weighted across events. Numerator excludes killed (`dn`);
  denominator is the untreated-net floor `(1 − rnm)` (bounded, ~0.76).
  *Previous formula `(sn+rnm)/sn` was buggy*: `/sn` blows up for insecticidal nets and inverted
  Pyr-ATN vs ATN exposure ordering; fixed by changing denominator to `(1 − rnm)`.
  No new parameters — existing bednet schedule matrices reused. If `match` returns NA (edge case:
  t0_atn not in bednet schedule), `contact_factor` falls back to 1 for that event. Tests in §12f.

## 11. Mali projection pipeline (dev/mali_projection_run.R, confirmed 2026-06-17)

- **Test harness:** `options(mali_test_mode = TRUE); source("dev/mali_projection_run.R")` loads all
  functions/data without launching the cluster. `dev/mali_single_region_test.R` uses this to run
  the highest-EIR region (Mopti, EIR ≈ 501) × 4 arms sequentially — use this to validate fixes
  before a full 36-job parallel run.
- **`site_parameters()` overrides clinical incidence rendering** with its own three age-group
  buckets: `[0, 1824]`, `[1825, 5474]`, `[5475, 36499]` (days). The `clinical_incidence_rendering_*`
  override in `render_overrides` is ignored. Sum the three columns for all-ages incidence.
- **Prevalence column off-by-one:** malariasimulation uses an exclusive upper bound in column names.
  `prevalence_rendering_max_ages = 10 * 365 = 3650` → column `n_detect_lm_730_3649` (not `_3650`).
  Same for `n_age_730_3649`. `run_one` uses these corrected names.
- **Mali species:** gambiae / arabiensis / funestus (3 species). Any parameter test with Mali must
  handle per-species vectors; single-species `get_parameters()` defaults keep those as scalars.
- **State sizes by arm in the Mali pipeline:** `none`/`cfp` use default `deltaq=1` (from
  `get_parameters()`; `atn_overrides` only sets `deltaq=deltaq_use` for ATN arms) → `n_states=27`;
  `atn`/`pyr_atn` use `deltaq=10` → `n_states=135`. Both sizes can fail at low mosquito density
  without the non-negativity floor + raised `a_tol`.
- **Bamako solver fix (implemented 2026-06-20):** `a_tol=0.1` (pipeline override) + C++ `nn`
  non-negativity floor in `create_eqs`. Without the floor, `a_tol=0.1` caused NaN cascades
  (negative-overshoot in dry-season trough → wrong-sign outflow → NaN). With both: confirmed
  the Erlang chain runs stably. See `dev/bamako_solver_diagnosis.md` for timings once available.
- **Validated run times (STALE — `n=1000`, `deltaq=10`, now using `n=10000`):** none ≈ 126 s,
  cfp ≈ 137 s, atn ≈ 229 s, pyr_atn ≈ 215 s. These were Mopti at human_pop=1000; current
  runs at 10000 are expected to take ~10× longer and are not yet benchmarked cleanly.
- **Net retention is site-sourced (changed 2026-06-18).** `retention_time` is no longer hardcoded;
  it reads `unique(site_obj$interventions$mean_retention)` (≈2014 d for MLI), with a top-of-script
  `retention_override <- NULL` hook for manual override. NB: the site value is far longer than the
  old hardcoded 588 d and materially changes the CD top-up coverage math (`cd_cov`).
- **Future net efficacy is resistance-projected per distribution year (changed 2026-06-18).**
  `build_future_schedule(arm, region)` (was `(arm, res)`) looks up projected pyrethroid resistance
  from `site_obj$vectors$pyrethroid_resistance` — a per-region, per-year table spanning **2000–2050**
  (NOT `interventions`, which only carries 2024 forward) — mapping each grid timestep to its calendar
  year (`start_year + grid %/% 365`) and calling `med_net(pars, res_year)`. So `dn0`/`rn`/`gamman`
  now vary across the future window in step with the rising resistance trend (verified Mopti:
  2025≈0.82 → 2031≈0.92). Past nets already read per-row site efficacy (unchanged).
- **`dn0_atn` override dropped for `pyr_atn` (2026-06-18).** Pyrethroid mortality for Pyr-ATN flows
  through the ITN-side `dn0`/`rn` in the net schedule (now resistance-projected); the ATN kernel's
  `dn0_atn` represents only the antimalarial's extra mortality (default 0) — avoids double-counting.
  See §3 "Repellency / pyrethroid-resistance coupling" for why repellency lives in `a`, not `delta_atn`.
- **Kernel invariant tests (2026-06-18):** `tests/testthat/test-atn-mosquito.R` §12d —
  six `compute_atn_kernels()` invariants. §12e (added 2026-06-18) — four `lambda_atn` auto-derive
  tests: auto-derive from `set_bednets(retention)`, explicit override respected, NULL fallback in
  kernel, logistic-retention warning + half-life mapping. Run after `devtools::load_all()` via
  `devtools::test_active_file()` (file must be the active RStudio tab) or select-all-and-run in editor.
  Resistance-sensitivity verification stays in `dev/mali_single_region_test.R` (needs full sim).
- **ODE solver error label (fixed 2026-06-18).** The "too much work" error in `src/solver.h` was
  previously hardcoded to "aquatic life stage model" regardless of which solver failed. The shared
  `Observer` now accepts a `model_name` string (passed as `"adult mosquito"` or `"aquatic mosquito
  larval"` from `create_adult_solver`/`create_aquatic_solver`). The error also now prints
  `n_states` to disambiguate (aquatic = 3; adult = 3 + (deltaq+1)*(2+spor_len)). Budget is
  per-day (observer.reset() each step). The adult solver carries the much larger state in this
  fork (27 non-ATN, 135 ATN-on states vs upstream 3) — so "too much work" on an ATN run is
  most likely adult, not aquatic.
