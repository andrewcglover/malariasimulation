# Task: finalise dev/atn_local_test.R against the live package API

I've added `dev/atn_local_test.R` — a stripped-down local sanity run for the ATN
model (2×2 EIR×resistance grid, four arms, median inputs). Please finalise it
against the actual package API by resolving the four `## VERIFY (Claude Code)`
markers:

1. The `get_parameters(overrides = list(...))` interface for the ATN params.
2. The rendering parameters that produce age-banded prevalence and clinical
   incidence (`prevalence_rendering_*`, `age_group_rendering_*`,
   `clinical_incidence_rendering_*` — or whatever the fork actually uses).
3. The exact rendered output column names — check `names()` on a real
   `run_simulation()` result and fix `pfpr2to10` / `clin_inc` / `timestep`
   accordingly.
4. The `set_bednets` and `run_simulation` argument names.

Also: change the `use_bompard` default in `get_parameters()` to TRUE and confirm
`use_eip_hill` already defaults to TRUE (both only affect ATN-on runs, so the
baseline tests are unaffected). Once done, the explicit `form_overrides` line in
the script can stay as a harmless safety net.

Constraints:
- Do **not** change the experimental design or the input-file paths
  (`dev/atn_params/…`, `dev/itn_params/…`).
- When it loads and runs cleanly, show me the two plots (prevalence-over-time and
  cases-averted) and flag anything that looks off.

Sanity expectation: the ATN / Pyr-ATN arms should outperform Pyr-CFP most at 90%
resistance (pyrethroid failing, antimalarial not). If ATN beats Pyr-CFP more at
high resistance than low, the mechanism is wired sensibly.
