# Bamako "too much work" / NaN ODE failure — diagnosis & fix

*Recorded 2026-06-20 · branch `atn-dev`*

---

## Symptom

Bamako (Mali admin-1 region with the lowest PfPR in the site file) fails with either:

```
Error: Solver error: too much work in the adult mosquito ODE solver.
```

or (with `a_tol = 0.1` and no non-negativity floor) completes but produces **NaN** outputs:

```
steps: 100000001, t: 5174.01, solver: adult mosquito, n_states: 27,
x[0]:7.27555e+10, x[1]:0.276537, x[2]:0.000719794,
x[3]:nan, x[4]:nan, ..., x[26]:nan
```

even with `ode_max_steps = 1e8`. The other eight Mali regions all complete cleanly.
`human_pop = 10000` (single smooth run).

---

## State sizes by arm

For the Mali pipeline, the `atn_overrides` block only sets `deltaq = deltaq_use` (10) for
the ATN arms; the `none` and `cfp` arms use the `get_parameters()` default `deltaq = 1`:

| Arm | `deltaq` | `spor_len` | `deltaqp1` | `n_states = 3 + deltaqp1*(2+spor_len)` |
|-----|---------|-----------|-----------|---------------------------------------|
| none, cfp | 1 | 10 | 2 | 27 |
| atn, pyr_atn | 10 | 10 | 11 | 135 |

The observed NaN cascade occurred on the `none` arm (27 states), confirming the issue is not
limited to the 135-state configuration.

---

## Why Bamako specifically

Bamako has the lowest EIR / PfPR in the site file. With `human_pop = 10000` and strong
Sahelian seasonality, the dry-season carrying capacity drops far enough that total adult
mosquito density approaches zero. Essentially: near-elimination in the mosquito population
for part of each year.

---

## Why the OLD 3-compartment model coped

Upstream malariasimulation's adult model (`Sm/Pm/Im`) integrated only **3 states** in the ODE.
Critically, the EIP was implemented as a **`std::deque<double> lagged_incubating`** of length
`tau ≈ 10` (days), **pushed and popped outside the ODE** in `adult_mosquito_model_update`. The
ODE derivative itself therefore contained no fast transit rate — its stiffest eigenvalue was
`μ ≈ 0.1 day⁻¹` (adult death rate). Near zero mosquitoes, the system decays gracefully as
`dS/dt ≈ −μ S` and the adaptive DOPRI5 stepper (boost `odeint`) sweeps an entire day in very few
internal steps.

---

## Why the NEW Erlang model fails there

The port moved the EIP into the ODE as an **Erlang chain** and added the **ATN exposure conveyor**,
introducing two fast rates that were previously outside the ODE:

| Rate | Value (defaults) | Source |
|------|-----------------|--------|
| `rho = spor_len / dem` | `10 / 10 = 1.0 day⁻¹` | EIP progression per Erlang stage |
| `kappa = deltaq / atn_window` | `1 / 10 = 0.1 day⁻¹` (none/cfp) or `1.0 day⁻¹` (atn) | ATN conveyor |

`rho` is ~10× faster than `μ`, making the Erlang system **stiffer** than the original.
Near zero, the fast rates continuously shuffle compartments among near-zero states.

### The two failure modes (both observed)

**1. Absolute-tolerance thrashing → "too much work".**
Boost's controlled DOPRI5 accepts a step when the local error satisfies
`err_i ≤ a_tol + r_tol · |x_i|` for all states. With `a_tol = 1e-4` (package default), the
solver must resolve sub-`1e-4`-mosquito quantities for all near-zero compartments while
`rho ≈ 1 day⁻¹` keeps shuffling them. Each rejected step halves `dt`; the budget is exhausted
within a single integration day.

**2. NaN cascade via negative-overshoot → NaN outputs (confirmed 2026-06-20).**
With `a_tol = 0.1`, the stepper can accept large steps that overshoot a near-zero compartment
to a small **negative** value. Unguarded negatives feed back into rate terms with the **wrong
sign** (outflow becomes inflow, creating positive feedback), which compounds until the
compartment blows up → NaN. The cascade then infects all adult states simultaneously (as seen
in the log: x[3..26] all NaN while x[0..2] aquatic states remain finite).

**`total_M` truncation is NOT the differentiator.** The
`static_cast<size_t>(total_M_d)` truncation exists identically in old and new models.

---

## Fix implemented (2026-06-20)

### Fix A — raise `a_tol` in the Mali pipeline

**IMPLEMENTED** in `dev/mali_projection_run.R` overrides: `a_tol = 0.1`.

This raises the absolute error floor from `1e-4` to `0.1` mosquitoes per compartment. The
accumulated error in `Ivtot` (e.g. 2 Iv compartments for the `none` arm) is at most ~0.2 —
well below the biological elimination floor of ~1 mosquito. Large compartments continue to be
controlled by `r_tol · |x_i|`, so accuracy is preserved where it matters.

**`a_tol = 0.1` not `1`:** the more aggressive `a_tol = 1` was considered but rejected because
accumulated Ivtot error (~O(10) for 11 compartments) could produce artefactual EIR spikes.
`a_tol = 0.01` was tried earlier but the run was killed prematurely at ~20 min — not confirmed
sufficient or insufficient.

**Properties:** no recompile, no signature change, one-line reverting edit.
Keep `r_tol = 1e-4` unchanged (irrelevant in the near-zero trough where `|x_i| ≈ 0`).

### Fix B — C++ non-negativity floor (IMPLEMENTED 2026-06-20)

**In `src/adult_mosquito_eqs.cpp`, `create_eqs` lambda:** all adult compartment reads are
wrapped through a `nn` (non-negative) helper lambda:

```cpp
auto nn = [&x](size_t i) -> double { return x[i] < 0.0 ? 0.0 : x[i]; };
```

All `x[sv_idx(...)]`, `x[ev_idx(...)]`, `x[iv_idx(...)]` and `x[AquaticState::P]` reads in
the derivative body use `nn(...)` instead of `x[...]` directly.

**Rationale:** Fix A (high `a_tol`) reduces solver work by allowing large steps, but those
large steps can overshoot near-zero to slightly negative, triggering Fix A's only failure mode
(NaN cascade). The floor prevents the cascade by returning `0` for any negative read —
breaking the positive-feedback loop before it can compound.

**Properties:**
- **Robustness, not speed.** The kink at zero may cost a few extra steps near a zero-crossing
  (DOPRI5 assumes smooth derivatives). This is a small, bounded overhead — not a speedup.
- Requires `devtools::load_all()` recompile; **no exported signature changes** (no
  `compileAttributes()` / RcppExports regen needed).
- Mass error from clamping is `O(a_tol)` per crossing — consistent with the already-accepted
  sub-`a_tol` erasure and biologically negligible (human `A`/`U` reservoir, not the
  short-lived mosquito pool, carries infection through the dry season).

### Biological rationale for accepting sub-`a_tol` erasure

Both fixes together numerically erase the fractional infectious-mosquito signal in the deep
trough. This is acceptable because:
- The **human asymptomatic/subpatent reservoir** (`A`/`U`, lasting weeks–months) is the
  realistic dry-season carrier — not the short-lived mosquito pool (`mu ≈ 0.1/day` ⇒ `Iv`
  decays e-fold per 10 d, effectively gone over a 3–4 month Sahelian dry season).
- The stochastic human IBM reseeds transmission whenever ≥ 1 human remains infected, independent
  of the mosquito `a_tol` setting.
- Genuine elimination requires all ~10 000 human individuals to clear in a trough — unlikely at
  Bamako's PfPR given the subpatent reservoir size.
- Reintroduction realism (metapopulation mixing) is explicitly **deferred** to a later workstream.

---

## Testing protocol (post-fix)

1. `devtools::load_all()` to recompile the C++ floor.
2. Run `dev/mali_single_region_test.R` with `test_region = "Bamako"` across all four arms.
3. Confirm all four arms complete without solver error and without NaN in outputs.
4. Sanity-check `dev/outputs/single_region_plots/`: dry-season near-zero EIR with wet-season
   recovery is expected and correct; persistent NaN or frozen-at-zero EIR would indicate a
   remaining issue.
5. Re-run the full 9-region × 4-arm grid to confirm no regressions in the eight passing regions.

---

## Escalation path (if Fix A + B insufficient)

### Fix C — higher-order explicit solver: RKF78

If Bamako is still slow post-floor (throughput, not stiffness, is the bottleneck), swap
`runge_kutta_dopri5` (4th/5th order, 6 fn evals/step) for `runge_kutta_fehlberg78` (7th/8th
order, 13 fn evals/step) in `src/solver.h`. Higher order means fewer steps to track sharp
seasonal transients. Requires `src/solver.h` edit + recompile.

### Fix D (last resort) — implicit solver (CVODE/Rosenbrock)

Only if stiffness ratio proves > ~100 (i.e. stability, not accuracy, is the constraint).
Both require the full 27×27 or 135×135 Jacobian and significant C++ work.

---

## Open questions

- Is sub-1-mosquito precision required for any planned output metric? (Prevalence, EIR, and
  incidence are aggregated over humans; sub-1-mosquito noise cannot propagate meaningfully.)
- Should `a_tol` be tuned globally or remain a pipeline override? Currently a pipeline-only
  override; the package default `1e-4` is correct for the stock 3-state model.
