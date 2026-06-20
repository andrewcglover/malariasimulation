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

Add `a_tol = 1` (or tune 0.1–1) to the `site_parameters()` overrides in
`dev/mali_projection_run.R:239-240` alongside `ode_max_steps`:

```r
overrides = c(..., list(ode_max_steps = 1e8, a_tol = 1))
```

**Rationale:** mosquito compartment counts are biologically ≥ 0 with elimination at ~1 mosquito.
Resolving sub-1-mosquito quantities to `a_tol = 1e-4` is physically meaningless and is the
proximate cause of the near-zero thrashing. With `a_tol = 1`, near-zero compartments are
allowed `O(1)` error (a mosquito or less), which is both biologically correct and eliminates the
solver death spiral. Large compartments continue to be controlled by `r_tol · |x_i|`, so
accuracy is preserved where it matters.

**Properties:** no recompile, no signature change, one-line reverting edit. Try `a_tol = 0.1`
first (conservative); escalate to `1` if still failing. Keep `r_tol = 1e-4` unchanged.

### Fix 2 (secondary, if Fix 1 is insufficient) — non-negativity floor in C++

Add a clamp in `src/adult_mosquito_eqs.cpp`, `create_eqs`, reading the compartment values
through `std::max(x[i], 0.0)` before computing derivatives. This prevents transiently negative
states from destabilising the stepper. Requires recompile + `Rcpp::compileAttributes()` + devtools
reload.

### Fix 3 (alternative, not recommended for production) — raise `ode_max_steps` further

Already at `1e8`; going to `1e9` would add ~10× wall time. Not a real fix — just postpones the
failure to a slightly lower population.

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
