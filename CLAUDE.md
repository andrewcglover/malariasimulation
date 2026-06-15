# CLAUDE.md — malariasimulation fork: add ATN (antimalarial-treated net) functionality

> Place this file in the **root of the malariasimulation fork** (`you/malariasimulation`, branch `atn-dev`).
> Claude Code reads it automatically each session. Keep it updated as decisions change.

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
   - **drop the `dt` factors** and the `if (x + dx < 0) 0` clamps — those are discrete-Euler hacks;
     the continuous adaptive solver doesn't need them.
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
- Continuous ODE ≠ discrete Euler: remove `dt` and negativity clamps; trust the solver tolerances.
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

**R-side touch-points (filenames still TO CONFIRM — find these next):**
- `get_parameters()` — add ATN params + no-ATN defaults.
- **mosquito equilibrium / init — `set_equilibrium()`**: the *human* equilibrium comes from the external
  `malariaEquilibrium` package (via `eq_params`) and is **unchanged** by ATN (it switches on later at
  `t0_atn`) — do **not** fork or modify `malariaEquilibrium`. The *mosquito* equilibrium is
  malariasimulation's own glue: `set_equilibrium()` derives the mosquito density + `init_foim` from that
  output, and the code building the `init` vector for `create_adult_solver` must now emit the enlarged
  `Sv/Ev/Iv` vector (whole equilibrium into baseline index 1 via `Ev_ratio = rho/(rho+mu)` &
  `Ev_norm_factor`, exposed rows zero). Because `malariaEquilibrium` — hence `init_foim` / total density —
  is unchanged, baseline init just *redistributes the same total*; but it reproduces stock only **to
  tolerance** (single-delay EIP vs Erlang; converges as `spor_len` grows), not bit-identically.
- the per-timestep process calling `adult_mosquito_model_update` — compute & pass the kernels.
- the infectious read-out feeding EIR/biting — swap the single `I` index for `sum(Iv)`.

**Compatibility strategy:** disaggregate internally, **sum only at the boundary** — `sum(Iv)` for the
EIR read-out and `sum(Sv)/sum(Ev)/sum(Iv)` for any `Sm/Pm/Im` outputs — and route `foim` specifically
into `Sv[1]` (exposed rows use `Lambda_i`, not baseline `foim`).

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
