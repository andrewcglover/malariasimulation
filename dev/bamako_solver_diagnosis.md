# Bamako "too much work" ODE failure — diagnosis

*Recorded 2026-06-20 · branch `atn-dev` · no code changes in this doc*

---

## Symptom

Bamako (Mali admin-1 region with the lowest PfPR in the site file) fails with:

```
Error: Solver error: too much work in the adult mosquito ODE solver.
Check for extreme changes in carrying capacity possibly as a result of
seasonality or set_carrying_capacity() parameterisation.
```

even with `ode_max_steps = 1e8` (raised in commit `9af1200`). The other eight
Mali regions all complete. `human_pop = 10000` (single smooth run).

---

## Why Bamako specifically

Bamako has the lowest EIR / PfPR in the site file. With `human_pop = 10000` and strong
Sahelian seasonality, the dry-season carry-capacity drops far enough that total adult
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
| `kappa = deltaq / atn_window` | `10 / 10 = 1.0 day⁻¹` | ATN exposure conveyor |

These are ~10× faster than `μ`, making the Erlang system **stiffer** than the original.

At the same time the state vector **expanded dramatically**:

| Arm | State size (`3 + deltaqp1*(2+spor_len)`) |
|-----|------------------------------------------|
| old (all arms) | 3 |
| new, none/cfp (baseline-only active) | `3 + 11*12 = 135` |
| new, atn/pyr_atn | same 135 (ATN-on activates all compartments) |

`src/adult_mosquito_eqs.cpp`, `create_eqs` loops over all 135 states every step.

### The two near-zero failure modes

**1. Absolute-tolerance thrashing.**  
Boost's controlled DOPRI5 accepts a step when the local error satisfies
`err_i ≤ a_tol + r_tol · |x_i|` for all states (see `src/solver.h`).
Both tolerances are `1e-4` by default (`R/parameters.R:547-548`).  
For the 100+ compartments that sit near zero during the dry season, the error floor is
`a_tol = 1e-4`. The solver must resolve sub-`1e-4`-mosquito quantities to this precision
while the fast rates `rho` and `kappa` (~1 day⁻¹) are continuously shuffling them.
Each rejected step halves the internal `dt`; the "too much work" budget (`ode_max_steps = 1e8`)
is exhausted within the one-day integration window.

**2. No negativity clamps on the continuous path.**  
The original discrete-Euler malsimple used `if (x + dx < 0) 0` guards. The port explicitly
dropped these for the continuous ODE (CLAUDE.md §8). Near zero, fast rates can transiently push
a compartment slightly negative; the adaptive stepper then shrinks `dt` aggressively to satisfy
the tolerance on what is effectively a stiff boundary at zero, compounding problem 1.

**`total_M` truncation is NOT the differentiator.** The line
`model.growth_model.total_M = static_cast<size_t>(total_M_d)` in `create_eqs` exists
identically in both the old and new models (it rounds the real-valued mosquito sum to an
integer for the larval model). This cannot explain why the new model fails; the differentiator
is the fast in-ODE rates × large near-zero state count × tight absolute tolerance.

---

## Recommended fixes (to apply after review)

### Fix 1 (primary) — raise adult ODE `a_tol` in the Mali pipeline

**IMPLEMENTED 2026-06-20** in `dev/mali_projection_run.R` overrides: `a_tol = 0.01`.

```r
overrides = c(..., list(ode_max_steps = 1e8, a_tol = 0.01))
```

**Rationale:** mosquito compartment counts are biologically ≥ 0 with elimination at ~1 mosquito.
Resolving sub-`1e-4`-mosquito quantities at the default `a_tol = 1e-4` is physically meaningless
for near-zero compartments. With `a_tol = 0.01`, the per-compartment absolute error floor rises
to 0.01 mosquitoes; the maximum accumulated error in `Ivtot` (11 Iv compartments) is ~0.1 —
well below the biological elimination floor of ~1 mosquito. Large compartments continue to be
controlled by `r_tol · |x_i|`, so accuracy is preserved where it matters.

`a_tol = 0.01` not `1`: the more aggressive `a_tol = 1` was considered but rejected because
accumulated Ivtot error (~O(10)) could produce artefactual EIR spikes during near-elimination.

**Properties:** no recompile, no signature change, one-line reverting edit.
Keep `r_tol = 1e-4` unchanged.

### Fix 2 (if Fix 1 insufficient) — higher-order explicit solver: RKF78

If `a_tol = 0.01` reduces but does not eliminate failures, the remaining cause is likely
**accuracy under the sharp dry-to-wet seasonal burst** (Bamako's carrying capacity rises steeply,
driving a rapid cascade: K spike → larval pool → pupae → `betaa` surge → all 10 Erlang stages
simultaneously). This is an accuracy-on-transient problem, not stiffness per se (Jacobian
eigenvalue ratio stays ~10–15).

Fix: swap `runge_kutta_dopri5` (4th/5th order, 6 fn evals/step) for
`runge_kutta_fehlberg78` (7th/8th order, 13 fn evals/step) in `src/solver.h`. Higher order means
the local error is O(h^7) vs O(h^5) — roughly 3–5× fewer steps to track the same transient.
Requires a small `src/solver.h` edit + recompile; no Jacobian or parameter changes.

*Do NOT try an implicit solver (Rosenbrock4/BDF) for this reason alone*: implicit methods
help when stability constrains step size (stiffness ratio > ~100); here the constraint is
accuracy, and implicit methods do not improve accuracy per se. They would also require the full
135×135 Jacobian (~135 extra fn evals per step) and C++ refactoring.

### Fix 3 (if Fix 2 insufficient) — non-negativity floor in C++

Add a clamp in `src/adult_mosquito_eqs.cpp`, `create_eqs`, reading near-zero compartments
through `std::max(x[i], 0.0)` before computing derivatives. Prevents transiently negative states
from destabilising the stepper on the recovery side of the trough. Requires recompile.

### Fix 4 (last resort) — implicit solver (CVODE/Rosenbrock)

If Fixes 1–3 all fail, the system is stiffer than the Jacobian analysis suggests (possibly
through indirect aquatic↔adult coupling effects). CVODE (Sundials) or boost's `rosenbrock4`
would then be appropriate. Both require the Jacobian (135×135, sparse-banded) and significant
C++ work. Flag for redesign if reached.

---

## Testing protocol (after applying Fix 1)

1. Run `dev/mali_single_region_test.R` configured for Bamako (not Mopti) across all four arms.
2. Confirm all four arms complete without error.
3. Spot-check that Bamako PfPR/EIR trajectories are physically sensible (zero-crossing and
   recovery in the dry season is expected and correct; persistent negatives would indicate a
   remaining issue).
4. Re-run the full 9-region × 4-arm grid to confirm no regressions in the previously passing
   eight regions.

---

## Open questions for the user

- Should `a_tol` apply globally (i.e. change the `get_parameters()` default from `1e-4`) or
  only in the Mali pipeline override? Given that stock malariasimulation uses 3 states and has
  no near-zero thrashing issue, the default is fine for upstream users — the override is safer.
- Is sub-1-mosquito precision required for any planned output metric? (Prevalence, EIR, incidence
  are all aggregated over the human population; mosquito compartment noise at the sub-1 level
  cannot propagate meaningfully.)
