# Task: make pre-infection blocking bulletproof — derive Lambda00 from B_max_post

## Cause (confirmed)

Pre-infection blocking is off because `Lambda00sf` defaults to `1` in the fork
(`R/parameters.R:490`), giving peak blocking `1 − 1 = 0`. The malariasimple
default is `0.0019` (which equals `1 − B_max_post` at the default
`B_max_post = 0.9981`). The malariasimple `get_parameters.R` is in
`dev/reference/get_parameters.R` for cross-reference.

This was also fragile in the original workflow: `B_max_post` is drawn per
posterior-draw, but `Lambda00sf` was left at its fixed default, so the two
**decoupled** for every draw except the median.

## Design decision: don't store `Lambda00sf` — derive it

Pre- and post-infection blocking share the same fitted max activity `b_max`:

- peak pre-infection blocking  = `1 − Lambda00sf`
- peak post-infection blocking = `B_max_post`  (= `b_max`)

so the identity `Lambda00sf = 1 − B_max_post` should **always** hold. Rather than
keep `Lambda00sf` as a free parameter that can fall out of sync (the bug above),
derive it at the **point of use** so it always reflects the current `B_max_post`
— median, override, or posterior draw, by any code path.

## Changes

1. **`R/biting_process.R`, `compute_atn_kernels()`** — replace the `Lambda00`
   computation (adapt to the exact local names):
   ```r
   # before
   Lambda00 <- foim * parameters$Lambda00sf
   # after
   Lambda00 <- foim * (1 - parameters$B_max_post)   # pre- & post- share b_max
   ```

2. **`R/parameters.R`** — remove the `Lambda00sf` default (line ~490); it's no
   longer a model parameter.

3. **Grep the whole repo for `Lambda00sf` and update every reference**, because
   `get_parameters()` errors on unknown parameters, so any leftover override will
   break:
   - `tests/testthat/test-atn-mosquito.R` (~line 127) currently sets
     `Lambda00sf = 0.7` — remove that override and re-derive the test's expected
     values from `B_max_post` (peak pre-blocking = `B_max_post`).
   - `dev/atn_local_test.R` — ensure `atn_overrides` does **not** set `Lambda00sf`
     (it already sets `B_max_post`).
   - any other script/doc.

4. **Align `p_atn`** in `dev/atn_local_test.R`: change `p_atn = 1.0` → `p_atn = 0.9`
   to match malariasimple's default.

5. Confirm nothing reads `Lambda00sf` on the **C++** side (it shouldn't —
   `compute_atn_kernels` builds `Lambda_i` R-side and passes the vector into the
   C++ update). If confirmed, this is an **R-only change**: no
   `Rcpp::compileAttributes()` and no exported-signature change — a plain
   `devtools::load_all()` is enough to pick it up.

## Verify (report findings; don't change the experimental design)

- Print `Lambda_i` across the exposure compartments in the pure-ATN run: it should
  now sit well **below** baseline `foim` near each distribution (it was equal to
  `foim` before).
- Re-run `dev/atn_local_test.R`: the pure-ATN arm should avert substantially more
  and reduce PfPR more, with the ATN advantage largest at 90% resistance.
- Confirm the baseline (no-ATN) regression tests still pass unchanged — deriving
  `Lambda00` doesn't affect ATN-off runs (where `delta_atn = 0`).
- Report the before/after `Lambda_i` and the updated plots.

## Note on an override hatch

Dropping `Lambda00sf` entirely is cleanest given pre/post should always be equal.
If you later need to decouple them for a sensitivity analysis, reintroduce
`Lambda00sf` as an explicit *optional* override that, when supplied, replaces
`1 − B_max_post` — but default to the derived value.
