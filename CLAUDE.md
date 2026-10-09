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
the antimalarial on an attempted bite. A mosquito exposed to an ATN (rate `xi`, below) moves
into exposure compartment 2. `Q_atn_t` = total ATN coverage across distribution events.

**Exposure rate — DESIGN DECISION (revised 2026-10-08; implements SI eq:xi).**
The rate passed to C++ (as the `av_da` argument) is computed in `simulate_bites`
(`R/biting_process.R`) as `xi = f/(1 − Z) · Q0 · delta_atn · p_contact`, i.e. SI `ξ = f_A Q_A`:

- **`f/(1 − Z)`** is the attempt rate `f_A` (SI eq:fA): feeding cycles per day × expected attempts
  per cycle `1/(1 − Z)`. `f` (`blood_meal_rate`) and `Z` (`average_p_repelled`) come from the feeding
  cycle already in `simulate_bites`; `Z` is population-level, so it includes historical ITNs and IRS.
- **`Q0`** is the mosquito anthropophagy `parameters$Q0[[s_i]]`, the share of *attempts* on humans.
  **Trap:** inside `compute_atn_kernels` the local `Q0` is the ATN coverage `Q0_atn`; that is why
  `xi` is built at the call site, not in the kernel.
- **`delta_atn`** (`p_atn * phi_bednets[[s]] * Q_t`, i.e. `p_A Φ_B U`) is the unconditional
  per-attempt share used in `xi`. Since 2026-10-08 it is **not** passed to C++ itself; the
  `delta_atn` C++ slot now carries ε (block below).
- **`p_contact`** (`compute_atn_kernels`) is the coverage-weighted mean over distribution events of
  `p_C = sn + ω·rn` (`sn = 1 − rn − dn`, `ω = omega_atn`, default 0.9, see §10), so
  `Q0 · delta_atn · p_contact` is SI eq:Q_A including its sum over events. Fed survivors contact the
  net; a repelled mosquito contacts it with probability ω whatever repelled it; pyrethroid-killed
  (`dn`) are excluded. Non-insecticidal ATN: `p_C = 1 − (1 − ω)·rnm` (0.976) at every net age. Fresh
  Pyr-ATN: `sn + ω·rn0`, below the ATN value (Pyr-ATN < ATN ordering), tending to it as the
  insecticide wanes. No `set_bednets` call, or no bednet row matching `t0_atn`: falls back to 1.
- ATN-off (`delta_atn = 0`) ⇒ `xi = 0`: **baseline untouched**. Only the exposure rate changed:
  `foim` (→ `Sv[0]`), EIR, `mu`/`f`, the aquatic model and total density are unaffected.
- Debug renders (`atn_debug = TRUE`): `dbg_p_contact_*` and `dbg_xi_*` (were `dbg_contact_factor_*`,
  `dbg_av_da_*`); read by `dev/segou_atn_debug.R`.

**History (replaced outright 2026-10-08, no switch; old behaviour = commit `5fdaa29`).** Was
`av_da = a · delta_atn · contact_factor`, `contact_factor = p_C/(1 − rnm)` (2026-06-19; the first
version divided by `sn`, which inverted the Pyr-ATN < ATN ordering). Two errors: `a = Q·f_R` counts
human *meals* (`Q`) rather than *attempts* (`Q0`), and `1/(1 − rnm)` stood in for the retry inflation
as if every attempt were on an ATN user. Ratio new/old `= (Q0/Q)(1 − rnm)/(1 − Z)`: −15.4% for pure
ATN at the illustrative config of `dev/reference/ATN_contact_factor_proposal.tex` (`U = 0.5`, ATN the
only net; the "~16%" quoted earlier was rounding). The proposal's `c = (χ/(d(1 − Z)))·p_contact`
gives the identical rate (`a·delta_atn·c`, the `Q` in `a` cancels `d = Q`). Every ATN-on result moves,
so the change rides on the single c24med re-run. Griffin 2010 SI Text S2 mapping verified 2026-07-22:
`f_R = 1/(δ1 + δ2)`, `δ1 = δ10/(1 − Z)`, `Q = 1 − (1 − Q0)/W`, `α = Q·f_R`; code `blood_meal_rate` =
`f_R`, `average_p_successful` = `W`, `average_p_repelled` = `Z`, `.human_blood_meal_rate` = `α`
(Griffin `.doc` eqns are embedded MathType — read prose via `antiword`).

**Dosed-infection share ε — DESIGN DECISION (2026-10-08; implements SI eq:varepsilon).**
ε is the probability that a successful feed on a human results in exposure. It is passed to C++ in
the `delta_atn` slot; every C++ use of `da` is the ε role (`(1 − da)·Λ_i` on each `Sv`, and
`da·Λ0_t·Svtot` into the first exposed latent stage). Built in `simulate_bites` as
`eps = delta_atn · s_feed / sum(.pi · prob_bitten_survives)`:
- **`s_feed`** (`compute_atn_kernels`) is the coverage-weighted mean over distribution events of
  `s_N = 1 − rn − dn`, from the same per-event `sn_e` as `p_contact`, with the same fallbacks (1).
- The denominator is `Σ_h π_h w_h`, which equals `(W − (1 − Q0))/Q0` by SI eq:W; it is computed
  directly (no cancellation, no division by `Q0`), and so includes historical ITNs and IRS.
- ATN-off ⇒ `eps = 0` exactly. Debug renders `dbg_s_feed_*`, `dbg_eps_*`.
- **History:** until 2026-10-08 (commit `0b0776d` and earlier) C++ received `delta_atn = p_A Φ_B U`,
  i.e. P(attempt is on an ATN user) rather than P(dosed | infecting feed). New/old `= s_feed/Σπw`:
  −15.4% for pure ATN at the illustrative config (`U = 0.5`, `s_N = 0.76`). That equals the ξ
  correction there because `W = 1 − Z` when nothing kills, so both reduce to `s_N/Σπw`. The drop is
  much larger for co-treated nets while the insecticide is fresh.

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
3. **Extra mortality** `dn_atn` — a `(1 - dn_atn)` survival factor on exposure-compartment *inflow*.
   **CORRECTED framing (2026-07-23):** this is NOT an antimalarial effect (antimalarials cause **no**
   excess mortality). It is the **pyrethroid** net mortality (`dn0`): set `dn0_atn = dn0` for an
   **AITN** (co-treated pyrethroid+antimalarial net), `0` for a pure ATN. **Held at 0 for all arms**
   to avoid double-counting — see the §10 `dn0_atn` note. The SI drops it from the drug-mechanism
   list (four → three).
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
- **`p_contact` sourcing (was `contact_factor`; denominator dropped 2026-10-08, see §3).**
  `compute_atn_kernels` reads `rn`/`rnm`/`dn0`/`gamman` for each ATN event via
  `match(t0_atn, parameters$bednet_timesteps)`. In the Mali pipeline every ATN event matches exactly
  one bednet row (same `timesteps` vector, `mali_projection_run.R:248-258`). Per event
  `p_C = sn_e + ω·rn_e`, `sn_e = 1 − rn_e − dn_e`, then coverage-weighted across events;
  no denominator (the retry inflation is `1/(1 − Z)` in `xi`). If `match` returns NA (edge case:
  t0_atn not in bednet schedule), `p_C` falls back to 1 for that event. `s_feed` (coverage-weighted
  `sn_e`, for ε) uses the same rows and the same fallbacks. Tests in §12f.
- **`omega_atn` = SI ω, P(net contact | repelled) (added 2026-10-09; replaced `chem_dose_atn`).**
  Default **0.9** in `get_parameters`, plausible range 0.8 to 1.0 (video tracking: Parker 2015,
  Gleave 2023, contact shares near-flat across untreated and treated nets). Applies to every
  repelled mosquito, barrier or insecticide, all net types and ages, so it moves **every** ATN arm,
  including pure `atn` (`p_C` 1 → 0.976). Fresh Pyr-ATN `p_C` more than doubles (Churcher 2024
  pyr-only at resistance 0.82: 0.303 → 0.717), since the old default said insecticide-repelled
  mosquitoes never contact the net. R-only. Pipeline: env var `ATN_OMEGA`; `OUT_FILE` gains
  `_omega{tag}` only when ω ≠ 0.9.
  **History:** `chem_dose_atn` (`f`, default 0, added 2026-06-23) gave `p_C = sn + rnm + f·(rn − rnm)`:
  every barrier-repelled mosquito contacted the net, only a fraction `f` of the insecticide-repelled
  did, so contact among repelled mosquitoes rose to 1 as the insecticide waned. Removed outright (no
  switch). The pipeline now **refuses** `ATN_CHEM_DOSE` if set, so the June driver scripts
  (`run_c24med_overnight.sh`, `run_c24med_rerun_DE.sh`), which set it, error instead of silently
  running the new model; they are kept unedited as the record of the June runs.
- **`atn_displace_t0` / `atn_displace_Q0` — non-drug displacement events (added 2026-06-23).**
  Supports mixed-delivery arms where a **non-drug net (e.g. Pyr-CFP mass campaign)** overwrites
  ATN holders in the IBM, causing `Q_atn_t` to collapse at each campaign. These parameters list
  the timesteps and coverages of such overwriting distributions; they enter `repl_factor` in
  `compute_atn_kernels` (`R/biting_process.R`) as `prod(1 − Q0_displace)` but contribute **nothing**
  to `Q_each`. Default `numeric(0)` (empty) → all existing arms bit-identical to pre-change
  behaviour. R-only, no C++ recompile. Used by the `pyr_cfp_mc_atn_cd` arm (see §12g tests,
  `dev/c24med_projection_run.R:build_params`). Key design points:
  - Only events **strictly later** than the drug event `i` and fired by `timestep` reduce its
    `repl_factor`; earlier displacement events do NOT reduce later drug events (correct chronology).
  - `cd_cov` / `cd_floor` math is **unchanged** — it defines the inter-campaign top-up target
    regardless of what product campaigns distribute. The sawtooth (build via CD → collapse at MC
    → rebuild) falls out naturally once `Q_atn_t` sees the displacement events.
  - `p_contact`, `lambda_atn` decay, `delta_atn`, Hill/Bompard kernels, C++ ODE: all unaffected.
- **`n_use_atn` render output (added 2026-06-23).** `compute_atn_kernels` returns `Q_t` (total ATN
  coverage at the current timestep). `simulate_bites` (`R/biting_process.R`) renders
  `n_use_atn = Q_t * human_population` once per timestep (gated on `s_i == 1L` since `Q_t` is
  species-independent) — the deterministic expected count of ATN net holders used by the mosquito
  ODE. Enables ATN-coverage time-series plots without instrumenting the IBM `net_time` variable.

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
- **`dn0_atn` = pyrethroid net mortality, held at 0 for all arms (CORRECTED 2026-07-23; supersedes
  the 2026-06-18 "dropped for pyr_atn" note).** `dn0_atn` is **NOT** the antimalarial's mortality
  (antimalarials cause none). It is the **pyrethroid** `dn0`, applied as a `(1 - dn_atn)` survival
  factor on the ATN exposure-compartment *inflow*; intended use is `dn0_atn = dn0` for an **AITN**
  (co-treated pyrethroid+antimalarial net), `0` for a pure ATN. Kept **0 for every arm** in
  production (`dev/c24med_projection_run.R` leaves it unset → 0). Confirmed correct to avoid
  **double-counting**: for `pyr_atn`, `set_bednets` gets the real pyrethroid `dn0`
  (`c24med:246,252,363-367`), which already elevates the death rate `mu` via `death_rate()`
  (`R/mosquito_biology.R:163-169`: `p1 = p1_0*W/(1 - Z*p1_0)`, `W` carries `sn = 1 - rn - dn`), so a
  `(1 - dn_atn)` at the exposure inflow would kill those deaths a second time.
  **Accepted limitation:** `mu` carries the pyrethroid mortality as a mean-field population-average,
  NOT concentrated on the net-contacting (antimalarial-exposed) mosquitoes — so exposed mosquitoes
  are not more likely to die, whereas for a true AITN the drug-contact and pyrethroid-kill are the
  same event (antimalarial effect for AITN arms is thus arguably slightly overestimated). Proper fix
  is structural (withhold net-contact mortality from `mu`, re-apply at the exposure transition,
  keeping `rn` driving `W`/`Z`) — deferred, limitation accepted 2026-10-08. **DONE (2026-10-08):**
  the four dev test scripts that previously set `dn0_atn = only$dn0` (`atn_local_test.R`,
  `atn_local_test_v2.R`, `atn_local_test_cd_v1.R`, `atn_convergence_check.R`) now set it to 0
  explicitly, so they no longer double-count if rerun.
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
- **IRS DDT `ms_gamma` sign error in MLI site file (found 2026-06-22).** The `site`-package
  MLI interventions table carries **`ms_gamma = +0.004` for every DDT IRS row (2000–2016)**
  — wrong sign. Correct behaviour requires `ms_gamma < 0` so `ms = spraying_decay(t, θ, γ)`
  decays to 0 as the insecticide wears off. With `+0.004`, `ms → 1` by ~5 yr, making
  `prob_spraying_repels → 1` (maximum repellency) for any recipient whose `spray_time` was set
  years ago — inflating Z and suppressing biting/ATN-exposure. Actellic rows (2017–2024) are
  correct (`ms_gamma = −0.009`). All 9 MLI regions are affected (153/225 rows); **Ségou** is
  worst (large DDT IRS 2008–2016, ~zero since 2016 → ~50-60% of pop with `rs → 1` in future
  window; Z_gambiae 0.69 vs ~0.19 elsewhere). **Mopti** unaffected (actellic IRS 2017–2024,
  correct sign). `spray_time` has no expiry mechanism (cf. bednets' `throw_away_net`), so
  stale values persist indefinitely.
  **Pipeline workaround:** `build_params()` in `dev/c24med_projection_run.R` applies
  `ms_ext$interventions$ms_gamma <- -abs(ms_ext$interventions$ms_gamma)` immediately after
  `expand_interventions()`, before `site::site_parameters()`. `-abs()` is idempotent (correct
  actellic rows are unaffected). Does NOT mutate `dev/site_files/without_split/MLI.rds`.
  **Provenance:** the IRS formula code (`spraying_decay`, `prob_spraying_repels`,
  `prob_bitten` spraying block) is 100% upstream (Giovanni Charles 2020–21); the fork never
  touched it. This is a site-data error, not a malariasimulation code bug.
  **TODO:** raise the DDT `ms_gamma` sign error with the `site`-package maintainer; consider
  an upstream warning in `set_spraying` for `ms_gamma > 0` (parked for now).
- **c24med pipeline env vars — `dev/c24med_projection_run.R` (confirmed 2026-06-24).** All are
  optional; defaults produce the standard multi-arm, site-retention output.

  | Env var | Effect | Default |
  |---|---|---|
  | `ANTIMAL_HL_YEARS` | Antimalarial half-life (years) → `hl{tag}` in `OUT_FILE` | `2.64` |
  | `ATN_OMEGA` | `omega_atn` (prob. a repelled mosquito contacts the net) → `_omega{tag}` suffix when ≠ 0.9 (see §10). `ATN_CHEM_DOSE` (removed 2026-10-09) now errors if set | `0.9` |
  | `SWEEP_ARMS` | Comma-separated arm subset (see §10) | all arms |
  | `SWEEP_CORES` | Parallel workers; pipeline scripts use `12` | `18` |
  | `OUT_SUFFIX` | Extra tag appended before `.rds` — keeps retention/add-on runs distinct | `""` |
  | `NET_RETENTION_DAYS` | Override site `mean_retention` (days); empty = site value (≈2014 d MLI) | site value |

  `OUT_FILE` pattern: `{iso}_c24med_projection_results_hl{hl_tag}{omega_suffix}{out_extra_suffix}.rds`
  (was `{chem_suffix}` before 2026-10-09). **Trap for the re-run:** a default run (ω = 0.9, no
  `OUT_SUFFIX`) gets the SAME filename as the June default outputs and overwrites them; set
  `OUT_SUFFIX` to keep the June results.
  Output metadata includes `retention_time` for traceability.
  Overnight driver: `dev/run_c24med_overnight.sh` (Jobs A–F sequential; June record, now errors on `ATN_CHEM_DOSE`).
  Targeted re-run: `dev/run_c24med_rerun_DE.sh` (Job D ATN-only + Job E; respects `SWEEP_CORES`).
- **Selectable multi-file plotter — `dev/c24med_projection_plots_select.R` (added 2026-06-24).**
  Per-series registry keyed on `(key, arm, f, hl, ret)`; non-ATN arms always pinned to the BASELINE
  job (Job C). Three selectable dimensions: `f` (`chem_dose_atn`, removed 2026-10-09; the `f`
  dimension is now stale, and the plotter is the user's to update), `hl`, `ret` (retention).
  Human-readable filename suffix (e.g. `kABCDEFGH_hl2p64_f0_ret1396`) + manifest CSV. RDS files
  cached per path; missing files warned + dropped (not error).
