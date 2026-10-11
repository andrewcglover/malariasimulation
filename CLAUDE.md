# CLAUDE.md — malariasimulation fork: ATN (antimalarial-treated net) functionality

> Lives in the **root of the malariasimulation fork** (`andrewcglover/malariasimulation`, branch
> `atn-dev`). Claude Code reads it automatically each session.
>
> **Maintenance (Claude Code):** keep this file current as you work. Whenever you confirm a source fact,
> locate a file / function / constructor signature, settle a design decision, or hit a non-obvious gotcha,
> update the relevant section concisely and replace any "TO CONFIRM" note with the finding — so future
> sessions inherit it. Keep edits tight; this is a working reference, not prose.

## 0. One-line goal

A fork of `mrc-ide/malariasimulation` with the **antimalarial-treated net (ATN)** mechanism in the
**deterministic adult-mosquito ODE**, kept an installable, backward-compatible R package that runs
on the Imperial DIDE cluster via `hipercow`.

**This repository holds the transmission model only.** Data, Stan fits, the 1000-row draws table and
the projection pipeline live in the analysis package **atnet** (`andrewcglover/atnet`, local
`C:\Users\ag4218\Local\GitHub\atnet`; projections in `analysis/07_run_projections.R`), which loads
this fork from source (`ATN_FORK_PATH`) or, on the cluster, installs it at a fixed tag.

## 1. Reference — READ THIS FIRST

The authority on the model is the **supplementary information (SI)** together with **this fork's
code**. The live SI is the user's file in OneDrive (currently
`Documents\ATN\Glover2026manuscript\09OCT26_ATN_Overleaf.tex`): read it, never edit it; copy it to a
scratch folder before working from it. §5 (combined model) gives the ODE system; §1 the exposure
rate ξ and the dosed-infection share ε; §4 the blocking kernels and the lab-to-field map; §6 the
uncertainty design.

The ATN mechanism was first ported (June 2026) from the user's discrete-time `malariasimple` model
(`malariasimple_aitn_deterministic_v3.R`). **That file is no longer authoritative**: the model has
changed substantially since (potency, exposure rate, dosed-infection share, contact probability).
Copies survive only in this repository's history (commit `967d2bd`, `dev/reference/`) and in the
user's `malariasimple_ATNs` repository.

## 2. What changes vs what stays

malariasimulation = **individual-based humans + deterministic compartmental ODE mosquitoes**
(aquatic `E/L/P` + adult compartments, integrated in C++ via the adult/aquatic solver classes).

**Only the adult mosquito infection block changed.** Everything else — aquatic stages, larval
seasonality, the entire human IBM, immunity, biting/EIR machinery — is upstream code.

## 3. The ATN mechanism

The 3-compartment adult model (`Sm/Pm/Im`) is replaced by a **2-D ATN-exposure structure**, in
**every** run (ATN on or off):

- `Sv[deltaqp1]`            susceptible, by ATN-exposure class
- `Ev[deltaqp1, spor_len]`  latent, ATN-exposure class × Erlang EIP stage
- `Iv[deltaqp1]`            infectious, by ATN-exposure class

Indices and rates:
- `deltaqp1 = deltaq + 1`. Class **1 = never exposed (baseline)**; **2..deltaqp1 = time since ATN
  exposure**. A conveyor at rate `kappa = deltaq / atn_window` (default `atn_window` 10 days, so
  `kappa` = 1/day at `deltaq` = 10) carries mosquitoes back towards baseline. Each class holds a
  full S, E^1..E^R, I set, so the adult state per species is `deltaqp1 × (2 + spor_len)`.
- `spor_len` (R, default 10) Erlang stages replace upstream's fixed EIP delay in **all** arms;
  baseline stage rate `rho = spor_len / dem` (1/day at `dem` 10).
- Midpoint convention: class `i` is read at `s_i = (i − 1.5) · atn_window/deltaq` days since exposure
  (1-based), i.e. SI `s_i = (i − ½)/κ`.

**Coverage is ATN-only by design.** `Q_atn_t` (`Q_t`) is built solely from the ATN distribution
events (`Q0_atn`/`t0_atn`, `n_atn` of them) in `compute_atn_kernels`, with exponential loss at
`lambda_atn` and random proportional replacement by later events (SI eq:cov_replace). It excludes
historical pyrethroid ITNs still in circulation — those suppress biting (via `net_time` in `a`) but
deliver no drug. Do **not** couple ATN coverage to the `W`/`Z` net-using population.

**Exposure rate ξ (SI eq:xi; implemented 2026-10-08, commit `0b0776d`).** Passed to C++ as the
`av_da` argument; computed in `simulate_bites` (`R/biting_process.R`) as
`xi = f/(1 − Z) · Q0 · delta_atn · p_contact`:
- **`f/(1 − Z)`** is the attempt rate `f_A` (SI eq:fA); `f` (`blood_meal_rate`) and `Z`
  (`average_p_repelled`) come from the feeding cycle, so `Z` includes historical ITNs and IRS.
- **`Q0`** is the anthropophagy `parameters$Q0[[s_i]]`. **Trap:** inside `compute_atn_kernels` the
  local `Q0` is the ATN coverage `Q0_atn`; that is why `xi` is built at the call site.
- **`delta_atn`** = `p_atn * phi_bednets[[s]] * Q_t` (SI `p_A Φ_B U`), unconditional per attempt.
- **`p_contact`** = coverage-weighted mean over events of `p_C = sn + ω·rn` (`sn = 1 − rn − dn`,
  `ω = omega_atn`). Fed survivors contact the net; a repelled mosquito contacts it with probability ω
  whatever repelled it; pyrethroid-killed (`dn`) are excluded. Pure ATN: `1 − (1 − ω)·rnm` at every
  age. No bednet row matching `t0_atn`: falls back to 1.
- ATN off (`delta_atn = 0`) ⇒ `xi = 0`, baseline untouched.
- History: until `5fdaa29` the rate was `a · delta_atn · p_C/(1 − rnm)` (meals rather than attempts,
  and `1/(1 − rnm)` in place of the retry inflation): −15.4% for pure ATN at `U = 0.5`. Griffin 2010
  Text S2 mapping: `f_R = 1/(δ1 + δ2)`, `δ1 = δ10/(1 − Z)`, `Q = 1 − (1 − Q0)/W`, `α = Q·f_R`; code
  `blood_meal_rate` = `f_R`, `average_p_successful` = `W`, `average_p_repelled` = `Z`,
  `.human_blood_meal_rate` = `α`.

**Dosed-infection share ε (SI eq:varepsilon; commit `de3f55e`).** The probability that a
successful feed on a human results in exposure; passed to C++ in the **`delta_atn` slot** (every
C++ use of `da` is ε: `(1 − da)·Λ_i` on each `Sv`, `da·Λ0_t·Svtot` into the first exposed latent
stage). `eps = delta_atn · s_feed / sum(.pi · prob_bitten_survives)`, where `s_feed` is the
coverage-weighted `s_N` (same rows and fallbacks as `p_contact`) and the denominator is `Σ_h π_h w_h`
(= `(W − (1 − Q0))/Q0`), so it includes historical ITNs and IRS. ATN off ⇒ `eps = 0`.

**Potency (SI eq:potency, eq:potency_mean; commits `6543691`, `5fdaa29`).** Event k's antimalarial
potency is `P_k = exp(−gamma_atn · age_k)`; `P_bar` is its coverage-weighted mean (0 when `Q_t = 0`).
One half-life governs all three mechanisms, applied on the **field (TBA) scale**, i.e. after the
Bompard map:
1. **Pre-infection blocking** `Lambda_i = Λ · (1 − P_bar · TBA_field(B_max_post · H_pre(s_i)))`
   (`Lambda_i[1] = Λ`); concurrent dose-and-infect `Lambda0_t = Λ · (1 − P_bar · TBA_field(B_max_post))`,
   the `s → 0` limit (SI eq:Lambda_0).
2. **EIP suppression** `rho_i = rho − (rho − rho0_t) · H_eip(s_i)`, with
   `rho0_t = rho − (rho − rho_frac · rho) · P_bar` (Hill form when `use_eip_hill`).
3. **Post-infection blocking** `B_post[j] = P_bar · TBA_field(B_max_post · H_post(t_post[j]))`,
   `t_post[j] = (j − 0.5) · dem/spor_len`: an infected mosquito re-exposed at stage j clears with
   that probability and moves to the first exposed class. Indexed by EIP stage only (agreed
   2026-09-21): a mosquito re-exposed after slowed development is treated as more recently infected
   than it was, which slightly overstates clearance.
`TBA_field` = Bompard map (`use_bompard`, `m_bompard` 1.57e-4, `k_bompard` 4.95e-6); kernel
parameters are lab-scale TRA values from the Stan fits (atnet draws table).

**`dn_atn` = `dn0_atn · P_bar` is inert.** It is a `(1 − dn_atn)` factor on exposure inflows meant
for an AITN's pyrethroid kill, NOT an antimalarial effect (antimalarials cause no excess mortality).
**Held at 0 for all arms** (closed 2026-10-08): the pyrethroid `dn0` already raises `mu` through
`death_rate()` (`p1 = p1_0·W/(1 − Z·p1_0)`), so applying it again would double-count. Accepted
limitation: `mu` carries that kill as a population average, not concentrated on dosed mosquitoes.

**Coupling:** human → mosquito FOI (`foim`) feeds `Sv[1]`; EIR reads `sum(Iv)`; the aquatic model
sees the sum of all adult compartments.

**Out-of-scope note:** the leMenach/Griffin death rate means a repellent-only net lowers `mu`
relative to no net, so mosquito density is slightly higher under ATN than no nets. Pre-existing
upstream behaviour; revisiting it is deferred.

## 4. Implementation (done, June 2026)

Route A, edit in place and keep it a package: the adult derivative and state vector rewritten to the
`Sv/Ev/Iv` system (`src/adult_mosquito_eqs.cpp/.h`), ATN parameters plumbed through
`get_parameters()` with no-ATN defaults, kernels computed R-side each step in `compute_atn_kernels`
(`R/biting_process.R`) and passed to C++, which only does the ODE arithmetic. Initialisation seeds
the baseline class at the no-ATN Erlang steady state with all exposed classes zero
(`initial_mosquito_counts()`, `R/mosquito_biology.R`). After C++ interface changes run
`Rcpp::compileAttributes()` then `devtools::load_all()`.

## 5. Verification

- **Baseline:** with ATN off the fork must reproduce stock malariasimulation EIR and prevalence to
  tolerance (Erlang EIP vs upstream's fixed delay, so not bit-identical).
  `tests/testthat/test-regression-sm-pm-im.R` compares against
  `tests/testthat/fixtures/reference_sm_pm_im.rds`, which is **git-ignored (local only)**; regenerate
  it with `dev/generate-reference.R` under the ORIGINAL package. `dev/check-equilibrium.R` checks the
  ATN-off equilibrium start.
- **ATN on:** checked against the SI equations, not against malariasimple.
  `tests/testthat/test-atn-mosquito.R` holds the kernel invariants (§12d), `lambda_atn` (§12e),
  `p_contact`/`s_feed` (§12f), displacement (§12g) and a mocked `simulate_bites` test of ξ; the
  C++ was read against SI §5 on 2026-10-08, and the SI §5.1 reductions hold as identities.
- Run after `devtools::load_all()`.

## 6. Environment, build & deployment

- **Local:** Windows + **Rtools** matching the R version; compile via `devtools::load_all()` /
  `devtools::install()`. Never install the fork over the user's stable malariasimulation.
- **Cluster:** Imperial **DIDE** cluster via **hipercow** (drivers `dide-windows`, `dide-linux`).
  atnet provisions the fork from GitHub at a fixed tag (`pkgdepends`:
  `andrewcglover/malariasimulation@<tag>`). C++ is portable, so either node type works.
  Re-provision after every C++ change.

## 7. Git workflow

- `origin` = `andrewcglover/malariasimulation`; `upstream` = `mrc-ide/malariasimulation`.
- `master` = stable line; `atn-dev` = active work. Merge `atn-dev → master` only when it compiles
  **and** passes the baseline test; **tag `master` at milestones** (e.g. `v0.1-atn`) so hipercow can
  provision exact versions.
- Sync upstream by merging `upstream/master → atn-dev` first (test), then `→ master`.

## 8. Watch-outs

- Keep backward compatibility: no-ATN defaults (`p_atn = 0`, `Q0_atn = 0`, `B_max_post = 0`,
  `rho_frac = 1`, `gamma_atn = 0`, `dn0_atn = 0`) must leave the model at its ATN-free dynamics.
- **Non-negativity floor (2026-06-20):** a `nn` lambda in `create_eqs` clamps adult compartment reads
  to `max(x[i], 0)`. Needed because a high `a_tol` (0.1 in atnet's pipeline) lets steps overshoot
  slightly negative, which otherwise cascades into NaN via wrong-sign outflows. See
  `dev/bamako_solver_diagnosis.md`.
- Don't rename the package — it would break the `malariasimulation::`-calling ecosystem.
- **Dimension bounds:** `deltaq ≥ 1`, `spor_len ≥ 1`, `n_atn ≥ 1`. For an ATN-off run use `n_atn = 1`
  with `Q0_atn = 0`; never `deltaq = 0`. Loops over `3:deltaqp1` / `2:spor_len` must be C-style and
  naturally empty at the minimum sizes.
- **Structural constants:** `deltaq`, `spor_len`, `n_atn` set array sizes at model construction, so
  a saved state can only be restored into a model of the same sizes.

## 9. Source facts (C++ and R entry points)

- `src/adult_mosquito_eqs.cpp` — derivative (`create_eqs`) + Rcpp exports
  `create_adult_mosquito_model`, `adult_mosquito_model_update`, `adult_mosquito_model_save_state`,
  `adult_mosquito_model_restore_state`, `create_adult_solver`. `src/adult_mosquito_eqs.h` — struct,
  index helpers, state size `3 + deltaqp1·(2 + spor_len)` (aquatic E/L/P at 0–2).
- Index layout: `Sv[q] → 3 + q`; `Ev[q,j] → 3 + deltaqp1 + q·spor_len + j`;
  `Iv[q] → 3 + deltaqp1·(1 + spor_len) + q` (0-based q, j).
- R side: `get_parameters()` (`R/parameters.R`) holds the ATN defaults; `set_equilibrium()` uses
  `malariaEquilibrium` for the human state only (do not modify `malariaEquilibrium`);
  `initial_mosquito_counts()` (`R/mosquito_biology.R`) emits the Erlang init vector;
  `R/compartmental.R` builds the solver (`kappa = deltaq/atn_window`) and `iv_block_indices()`;
  `simulate_bites` (`R/biting_process.R`) calls `compute_atn_kernels` and
  `adult_mosquito_model_update(model, mu, foim, xi, eps, dn_atn, Lambda0_t, Lambda_i, rho_i, B_post, f)`
  each timestep; EIR reads `sum(Iv)`.
- Names: `dem` = EIP duration (days), `mum` = adult death rate, `foim`/`init_foim` = FOI on mosquitoes.

## 10. Naming, collisions & parameter conventions

- **Only one shared namespace matters:** the `get_parameters()` list. R lists silently tolerate
  duplicate names, so check every new key against the existing list. Suffix ATN keys with `_atn`;
  keep C++ internals short. `Q0` (anthropophagy) vs `Q0_atn` (ATN coverage) — distinct.
- **ATN distribution events:** one code path for one or many rounds (`n_atn`, vectors `t0_atn`,
  `Q0_atn`); replacement mixing assumes `t0_atn` chronological and distinct.
- **`use_bompard` and `use_eip_hill` default TRUE.** Only affect ATN-on runs.
- **`Lambda00sf` was removed** (2026-06-16); passing it errors.
- **`p_atn`** (default 0) must be positive when distributing ATNs, or the mechanism is silently off.
- **`set_bednets` requires `rnm < rn` (strict).** Non-insecticidal ATN: `rnm = rn − 1e-9`.
- **`phi_bednets` is per-species:** use `parameters$phi_bednets[[species]]`; a bare vector gives
  `Expecting a single value: [extent=N]` from Rcpp. Call sites pass the species index.
- **`lambda_atn` (default NULL = auto-derive).** `set_bednets` sets `lambda_atn = 1/bednet_retention`
  (logistic retention: warning + `1/bednet_logistic_half_life`), the mean-field analogue of the IBM's
  Exponential(mean = retention) net loss (`log_uniform`, `R/utils.R`). An explicit value is respected;
  no `set_bednets` call ⇒ 0. **Trap:** it is derived on the FIRST `set_bednets` call only, and
  `site::site_parameters()` already calls it with the site file's retention. A pipeline that calls
  `set_bednets` again with another retention must set `lambda_atn` itself (atnet does, since
  `d26d824`).
- **`p_contact` / `s_feed` sourcing:** each ATN event's `rn`/`rnm`/`dn0`/`gamman` come from
  `match(t0_atn, parameters$bednet_timesteps)`; in atnet every ATN event has its own bednet row.
- **`omega_atn`** = SI ω, P(net contact | repelled), default 0.9; applies to every repelled mosquito,
  all net types and ages (replaced `chem_dose_atn` 2026-10-09). In atnet ω is drawn per row
  (Beta(41.8, 4.65)).
- **`atn_displace_t0` / `atn_displace_Q0`:** non-drug distributions that overwrite ATN holders;
  they reduce `repl_factor` but add nothing to coverage. Default empty. Used only by the mixed arm
  `pyr_cfp_mc_atn_cd`, which is not in the paper.
- **`n_use_atn` render:** `Q_t × human_population`, the expected number of ATN users (once per
  timestep).
- **Debug renders** (`atn_debug = TRUE`): `dbg_a_`, `dbg_W_`, `dbg_Z_`, `dbg_f_`, `dbg_delta_atn_`,
  `dbg_p_contact_`, `dbg_xi_`, `dbg_s_feed_`, `dbg_eps_` per species.

## 11. Running the fork in projections (pipeline itself is in atnet)

Facts about the model that matter when running it (all confirmed in the Mali work):
- **`site_parameters()` sets clinical incidence rendering** to its own three age groups `[0, 1824]`,
  `[1825, 5474]`, `[5475, 36499]` days, ignoring `clinical_incidence_rendering_*`; sum the three.
- **Rendered column names use an exclusive upper bound:** `prevalence_rendering_max_ages = 3650` gives
  `n_detect_lm_730_3649` and `n_age_730_3649`.
- **Mali has three species** (gambiae, arabiensis, funestus): per-species vectors everywhere.
- **State sizes:** 27 per species at `deltaq` = 1 (non-ATN default), 135 at `deltaq` = 10 (ATN arms).
  Both can fail at low mosquito density without the `nn` floor + `a_tol = 0.1` (§8).
- **Solver error label:** the "too much work" error names the model (`adult mosquito` / `aquatic
  mosquito larval`) and prints `n_states`; on ATN runs it is most likely the adult solver.
- **MLI site file DDT IRS `ms_gamma` sign error:** every DDT IRS row (2000–2016) has
  `ms_gamma = +0.004` (should be negative), so repellency rises to 1 instead of decaying; actellic
  rows (2017+) are correct. Worst in Ségou. atnet's `build_params()` applies `-abs(ms_gamma)` after
  `expand_interventions()`. Upstream IRS code is untouched; it is a site-data error. TODO (parked):
  raise it with the `site` maintainer.
