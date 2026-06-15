# ATN port — step-by-step implementation plan

Route A: edit in place, keep the package installable and backward-compatible.

> **Atomicity note (updated):** Step 1 is independently completable (R only, no compile needed).
> **Steps 2–9 form a single atomic "compile-and-wire" block.** The C++ interface change
> (Steps 2–3) immediately breaks the R call sites; Steps 4–9 fix them. The package will
> not compile at any intermediate point between Steps 2 and 9. Do not attempt a test run
> until Step 9 is complete and `Rcpp::compileAttributes()` + `devtools::load_all()` succeed.
> **First test checkpoint: after Step 9** (behavioural smoke test run by the user).
> Steps 10–12 (rendering, final compile, verification tests) follow the smoke test.

---

## Step 1 — Add ATN parameters to `get_parameters()` (`R/parameters.R:341`)

Add the entries below inside `get_parameters()`, grouped after the existing
bednets block. Every value is a no-ATN default so nothing changes for callers
that don't set them.

**Structural constants** (fix array sizes at model-construction time):

```r
deltaq   = 1L,     # ATN-exposed compartments; deltaqp1 = deltaq + 1
spor_len = 10L,    # Erlang stages approximating EIP
n_atn    = 1L,     # number of ATN distribution events
```

**Per-event vectors** (length `n_atn`):

```r
t0_atn  = 0,       # day of each distribution event
Q0_atn  = 0,       # initial coverage of each event
```

**Scalar decay / mixing parameters**:

```r
lambda_atn  = 0,   # ATN coverage retention decay rate
gamma_atn   = 0,   # drug-effect potency decay rate
p_atn       = 0,   # probability antimalarial is present on bite attempt
```

**Drug-effect baselines**:

```r
Lambda00sf  = 1,   # scale factor: Lambda00 = Lambda * Lambda00sf
rho_frac    = 1,   # scale factor: rho00 = rho_frac * rho
dn0_atn     = 0,   # peak extra mortality
```

**Hill kernel parameters** (pre-infection blocking, EIP suppression,
post-infection blocking):

```r
s_half_pre  = 1,  nH_pre  = 1,
s_half_eip  = 1,  nH_eip  = 1,
B_max_post  = 0,  s_half_post = 1,  nH_post = 1,
```

**Bompard transform switches** (applied to lab→field conversion):

```r
use_bompard  = FALSE,
use_eip_hill = TRUE,
m_bompard    = 1.57e-4,
k_bompard    = 4.95e-6,
```

No existing key conflicts: `Q0` (anthropophagy) is distinct from `Q0_atn`.

---

## Step 2 — Rewrite `src/adult_mosquito_eqs.h`

### 2a. Replace the `AdultState` enum with index helper functions

Remove:
```cpp
enum class AdultState : size_t {S = 3, E = 4, I = 5};
```

Add (inline, after the aquatic include):
```cpp
// Index helpers — call after constructing the model (deltaqp1, spor_len known)
inline size_t sv_idx(size_t q)            { return 3 + q; }
inline size_t ev_idx(size_t q, size_t j,
                     size_t deltaqp1)     { return 3 + deltaqp1 + q * spor_len + j; }
inline size_t iv_idx(size_t q,
                     size_t deltaqp1,
                     size_t spor_len)     { return 3 + deltaqp1 * (1 + spor_len) + q; }
inline size_t state_size(size_t deltaqp1,
                         size_t spor_len) { return 3 + deltaqp1 * (2 + spor_len); }
```

These replace the `get_idx(AdultState::X)` pattern throughout the `.cpp`.

### 2b. Rewrite the `AdultMosquitoModel` struct

Remove: `std::deque<double> lagged_incubating`, `const double tau`.

Add:
```cpp
struct AdultMosquitoModel {
    AquaticMosquitoModel growth_model;
    // structural constants (fixed at construction)
    size_t deltaq, deltaqp1, spor_len;
    double kappa;   // = 1.0 / deltaq
    double rho;     // = spor_len / dem  (baseline EIP rate)
    // per-step scalars (updated by adult_mosquito_model_update)
    double mu, foim, delta_atn, dn_atn;
    // per-step vectors (length deltaqp1 or spor_len)
    std::vector<double> Lambda_i;   // length deltaqp1
    std::vector<double> rho_i;      // length deltaqp1
    std::vector<double> B_post;     // length spor_len
    AdultMosquitoModel(AquaticMosquitoModel, size_t, size_t, double, double, double);
};
```

---

## Step 3 — Rewrite `src/adult_mosquito_eqs.cpp`

### 3a. Constructor

```cpp
AdultMosquitoModel::AdultMosquitoModel(
    AquaticMosquitoModel growth_model,
    size_t deltaq,
    size_t spor_len,
    double mu,
    double dem,         // parameters$dem = delayMos
    double foim
) : growth_model(growth_model),
    deltaq(deltaq), deltaqp1(deltaq + 1), spor_len(spor_len),
    kappa(1.0 / deltaq), rho(spor_len / dem),
    mu(mu), foim(foim),
    delta_atn(0.0), dn_atn(0.0),
    Lambda_i(deltaq + 1, foim),   // default: no ATN → Lambda_i[i] = foim
    rho_i(deltaq + 1, spor_len / dem),
    B_post(spor_len, 0.0)
{}
```

### 3b. `total_M` update in `create_eqs`

Replace the current three-term sum with:
```cpp
double total_M = 0.0;
for (size_t q = 0; q < model.deltaqp1; ++q) {
    total_M += x[sv_idx(q)];
    for (size_t j = 0; j < model.spor_len; ++j)
        total_M += x[ev_idx(q, j, model.deltaqp1)];
    total_M += x[iv_idx(q, model.deltaqp1, model.spor_len)];
}
model.growth_model.total_M = total_M;
```

### 3c. ODE derivative — transcribe from v3 lines 537–603, **without `dt` factors and without negativity clamps**

The continuous solver handles small steps; clamping `if (x + dx < 0) 0` is
wrong in a continuous ODE.

**dSv:**
```
dSv[0] = betaa + kappa*Sv[deltaq] - (av*delta_atn + (1-delta_atn)*Lambda_i[0] + mu)*Sv[0]
dSv[1] = delta_atn*(av*(1-dn_atn) - Lambda0_t)*Svtot
          + av*delta_atn*(1-dn_atn)*sum_j(B_post[j]*Ecol[j])
          - (av*delta_atn + (1-delta_atn)*Lambda_i[1] + kappa + mu)*Sv[1]
dSv[i] = kappa*Sv[i-1] - (av*delta_atn + (1-delta_atn)*Lambda_i[i] + kappa + mu)*Sv[i]
          for i = 2..deltaq
```
(`betaa = 0.5 * PL / dPL` comes from the aquatic coupling already in `growth_eqs`.)

**dEv:**
```
dEv[0,0] = kappa*Ev[deltaq,0] + (1-delta_atn)*Lambda_i[0]*Sv[0]
            - (av*delta_atn + rho + mu)*Ev[0,0]

dEv[1,0] = (1-B_post[0])*av*delta_atn*(1-dn_atn)*Ecol[0]
            + delta_atn*(1-dn_atn)*Lambda0_t*Svtot
            + (1-delta_atn)*Lambda_i[1]*Sv[1]
            - (av*delta_atn + kappa + rho_i[1] + mu)*Ev[1,0]

dEv[i,0] = kappa*Ev[i-1,0] + (1-delta_atn)*Lambda_i[i]*Sv[i]
            - (av*delta_atn + kappa + rho_i[i] + mu)*Ev[i,0]
            for i = 2..deltaq

dEv[0,j] = kappa*Ev[deltaq,j] + rho*Ev[0,j-1]
            - (av*delta_atn + rho + mu)*Ev[0,j]
            for j = 1..spor_len-1

dEv[1,j] = (1-B_post[j])*av*delta_atn*(1-dn_atn)*Ecol[j] + rho_i[1]*Ev[1,j-1]
            - (av*delta_atn + kappa + rho_i[1] + mu)*Ev[1,j]
            for j = 1..spor_len-1

dEv[i,j] = kappa*Ev[i-1,j] + rho_i[i]*Ev[i,j-1]
            - (av*delta_atn + kappa + rho_i[i] + mu)*Ev[i,j]
            for i = 2..deltaq, j = 1..spor_len-1
```

**dIv:**
```
dIv[0] = kappa*Iv[deltaq] + rho*Ev[0,spor_len-1] - (av*delta_atn + mu)*Iv[0]
dIv[1] = av*delta_atn*(1-dn_atn)*Ivtot + rho_i[1]*Ev[1,spor_len-1]
          - (av*delta_atn + kappa + mu)*Iv[1]
dIv[i] = kappa*Iv[i-1] + rho_i[i]*Ev[i,spor_len-1]
          - (av*delta_atn + kappa + mu)*Iv[i]
          for i = 2..deltaq
```

Note on indexing convention: v3 uses 1-based Erlang (`j=1..spor_len`); C++
uses 0-based (`j=0..spor_len-1`). `rho` applies at `j=0` (first Ev stage),
`rho_i` at exposed compartments. Baseline row uses `rho` throughout; exposed
rows use `rho_i[i]`.

### 3d. Rewrite `create_adult_mosquito_model`

New signature adds `deltaq`, `spor_len`, `dem`; removes `tau`, `susceptible`,
`incubating`:
```cpp
//[[Rcpp::export]]
Rcpp::XPtr<AdultMosquitoModel> create_adult_mosquito_model(
    Rcpp::XPtr<AquaticMosquitoModel> growth_model,
    double mu,
    size_t deltaq,
    size_t spor_len,
    double dem,
    double foim
)
```

### 3e. Rewrite `adult_mosquito_model_update`

Remove deque push/pop. Accept new per-step ATN vectors:
```cpp
//[[Rcpp::export]]
void adult_mosquito_model_update(
    Rcpp::XPtr<AdultMosquitoModel> model,
    double mu,
    double foim,
    double delta_atn,
    double dn_atn,
    std::vector<double> Lambda_i,
    std::vector<double> rho_i,
    std::vector<double> B_post,
    double f
)
```
Body: assign all fields; update `growth_model.f` and `growth_model.mum`.

### 3f. Simplify `save_state` / `restore_state`

No deque to checkpoint — the ODE solver already owns the state vector.
These functions can become no-ops or be removed (check whether the event-
system save/restore path calls them before removing).

### 3g. `create_adult_solver`

Body is unchanged; `init` is now longer (see Step 6) but the function just
passes it through to the `Solver` constructor.

After all C++ edits: run `Rcpp::compileAttributes()` to regenerate
`R/RcppExports.R`, then `devtools::load_all()` to confirm it compiles.

---

## Step 4 — Expand `ADULT_ODE_INDICES` and add Iv-block helpers (`R/compartmental.R:1–2`)

`ADULT_ODE_INDICES` is currently a hardcoded named integer vector. Replace it
with a function that takes `deltaq` and `spor_len` and returns a named vector
covering all compartments, **plus** a separate helper that returns the Iv-block
range for the EIR sum.

```r
make_adult_ode_indices <- function(deltaq, spor_len) {
  deltaqp1 <- deltaq + 1L
  base <- 3L  # aquatic occupies positions 1-3
  sv_idx <- base + seq_len(deltaqp1)                              # Sv[1..deltaqp1]
  ev_idx <- base + deltaqp1 + seq_len(deltaqp1 * spor_len)       # Ev (row-major)
  iv_idx <- base + deltaqp1 * (1L + spor_len) + seq_len(deltaqp1)# Iv[1..deltaqp1]
  c(
    setNames(sv_idx, paste0("Sv", seq_len(deltaqp1))),
    setNames(ev_idx, paste0("Ev", seq_len(deltaqp1 * spor_len))),
    setNames(iv_idx, paste0("Iv", seq_len(deltaqp1)))
  )
}

iv_block_indices <- function(deltaq, spor_len) {
  deltaqp1 <- deltaq + 1L
  3L + deltaqp1 * (1L + spor_len) + seq_len(deltaqp1)
}
```

Keep `ADULT_ODE_INDICES <- c(Sm = 4, Pm = 5, Im = 6)` as a legacy alias
returning the equivalent subset when `deltaq = 1, spor_len = 1` — or retire it
and update every call site. Call sites to update:

| Location | Current use | ATN replacement |
|---|---|---|
| `R/compartmental.R:51` | `[ADULT_ODE_INDICES['Sm']]` — initial susceptible | `[sv_idx(1)]` — baseline Sv[1] |
| `R/compartmental.R:99` | `c(ODE_INDICES, ADULT_ODE_INDICES)` — rendering | `c(ODE_INDICES, make_adult_ode_indices(...))` |
| `R/biting_process.R:208` | `solver_states[[ADULT_ODE_INDICES['Sm']]]` — foim update | `solver_states[[sv_idx(1)]]` |
| `R/biting_process.R:263` | `solver_states[[ADULT_ODE_INDICES['Im']]]` — EIR read-out | see Step 7 |

---

## Step 5 — Compute per-step ATN kernels in R

Add a new function (e.g. `compute_atn_kernels(t, parameters)`) called once per
timestep before `adult_mosquito_model_update`. This is the direct analogue of
the ITN net-decay computation in the bednets process.

**`delta_atn` and coverage-weighted effect terms** (v3 lines 773–864 → R):

```r
compute_atn_kernels <- function(t, parameters, foim) {
  n  <- parameters$n_atn
  t0 <- parameters$t0_atn
  Q0 <- parameters$Q0_atn

  # Per-event coverage with random proportional replacement
  repl_factor <- vapply(seq_len(n), function(i) {
    later <- which(t0 > t0[i] & t0 <= t)
    prod(1 - Q0[later])
  }, numeric(1))
  Q_each <- ifelse(t < t0, 0,
              Q0 * exp(-parameters$lambda_atn * (t - t0)) * repl_factor)
  Q_t <- sum(Q_each)

  # Per-event drug-effect decay
  rho    <- parameters$spor_len / parameters$dem
  Lambda <- foim
  Lambda00 <- Lambda * parameters$Lambda00sf
  rho00    <- parameters$rho_frac * rho
  age <- pmax(t - t0, 0)
  Lambda0_each <- ifelse(t < t0, Lambda,
                   Lambda - (Lambda - Lambda00) * exp(-parameters$gamma_atn * age))
  rho0_each    <- ifelse(t < t0, rho,
                   rho    - (rho    - rho00)    * exp(-parameters$gamma_atn * age))
  dn_each      <- ifelse(t < t0, 0,
                   parameters$dn0_atn * exp(-parameters$gamma_atn * age))

  # Coverage-weighted averages
  Lambda0_t <- if (Q_t > 0) sum(Q_each * Lambda0_each) / Q_t else Lambda
  rho0_t    <- if (Q_t > 0) sum(Q_each * rho0_each)    / Q_t else rho
  dn_atn    <- if (Q_t > 0) sum(Q_each * dn_each)      / Q_t else 0

  delta_atn <- parameters$p_atn * parameters$phi_bednets * Q_t

  deltaqp1 <- parameters$deltaq + 1L

  # Bompard helper
  bompard <- function(b_lab) {
    m <- parameters$m_bompard; k <- parameters$k_bompard
    a <- (k / (k + m))^k
    b <- (k / (k + m * (1 - b_lab)))^k
    (b - a) / (1 - a)
  }

  # Lambda_i (length deltaqp1)
  s <- seq_len(deltaqp1) - 0.5
  b_lab_pre <- (1 - Lambda0_t / Lambda) *
    (parameters$s_half_pre^parameters$nH_pre /
     (parameters$s_half_pre^parameters$nH_pre + s^parameters$nH_pre))
  b_field_pre <- if (parameters$use_bompard) bompard(b_lab_pre) else b_lab_pre
  Lambda_i <- Lambda * (1 - b_field_pre)

  # rho_i (length deltaqp1)
  rho_i <- if (parameters$use_eip_hill) {
    rho - (rho - rho0_t) *
      (parameters$s_half_eip^parameters$nH_eip /
       (parameters$s_half_eip^parameters$nH_eip + s^parameters$nH_eip))
  } else {
    rho - (rho - rho0_t) * exp(-parameters$zeta * s)
  }

  # B_post (length spor_len)
  t_post <- (seq_len(parameters$spor_len) - 0.5) * parameters$dem / parameters$spor_len
  b_lab_post <- parameters$B_max_post *
    (parameters$s_half_post^parameters$nH_post /
     (parameters$s_half_post^parameters$nH_post + t_post^parameters$nH_post))
  B_post <- if (parameters$use_bompard) bompard(b_lab_post) else b_lab_post

  list(delta_atn = delta_atn, dn_atn = dn_atn,
       Lambda_i = Lambda_i, rho_i = rho_i, B_post = B_post)
}
```

---

## Step 6 — Enlarge `initial_mosquito_counts()` (`R/mosquito_biology.R:10`)

Replace the current 6-element return with the full `3 + deltaqp1*(2+spor_len)`
vector. Keep the aquatic block (`n_E`, `n_L`, `n_P`) unchanged; replace the
adult block.

```r
initial_mosquito_counts <- function(parameters, species, foim, m) {
  # ... existing aquatic n_E, n_L, n_P calculations unchanged ...

  mum      <- parameters$mum[[species]]
  deltaq   <- parameters$deltaq
  deltaqp1 <- deltaq + 1L
  spor_len <- parameters$spor_len
  rho      <- spor_len / parameters$dem

  n_Sm <- m * mum / (foim + mum)

  # Erlang redistribution of n_Pm into Ev[1, 1..spor_len]
  n_Pm <- m * foim / (foim + mum) * (1 - exp(-mum * parameters$dem))
  Ev_ratio      <- rho / (rho + mum)
  Ev_norm_factor <- if (Ev_ratio == 1) spor_len
                    else (1 - Ev_ratio^spor_len) / (1 - Ev_ratio)
  E1_eq <- n_Pm / Ev_norm_factor
  Ev_baseline <- E1_eq * Ev_ratio^(seq_len(spor_len) - 1)  # Ev[1, j], j=1..spor_len

  n_Im <- m * foim / (foim + mum) * exp(-mum * parameters$dem)

  c(
    n_E, n_L, n_P,                           # aquatic (unchanged)
    n_Sm, rep(0, deltaq),                     # Sv[1..deltaqp1]: baseline=n_Sm, exposed=0
    Ev_baseline,                              # Ev[1, 1..spor_len]
    rep(0, deltaq * spor_len),               # Ev[2..deltaqp1, *] = 0
    n_Im, rep(0, deltaq)                      # Iv[1..deltaqp1]: baseline=n_Im, exposed=0
  )
}
```

Layout in the returned vector (1-based R indices):

| Slice | Contents | Length |
|---|---|---|
| `1:3` | n_E, n_L, n_P | 3 |
| `4:(3+deltaqp1)` | Sv[1..deltaqp1] | deltaqp1 |
| `(4+deltaqp1):(3+deltaqp1+deltaqp1*spor_len)` | Ev row-major | deltaqp1 × spor_len |
| `(4+deltaqp1*(1+spor_len)):(3+deltaqp1*(2+spor_len))` | Iv[1..deltaqp1] | deltaqp1 |

---

## Step 7 — Update `parameterise_mosquito_models()` (`R/compartmental.R:4`)

Pass `deltaq`, `spor_len`, `dem` to `create_adult_mosquito_model`:

```r
create_adult_mosquito_model(
  growth_model,
  parameters$mum[[i]],
  parameters$deltaq,
  parameters$spor_len,
  parameters$dem,
  parameters$init_foim
)
```

The `susceptible * foim` initialisation of the deque is gone; the enlarged
`init` vector (from `initial_mosquito_counts`) now carries the full equilibrium.

---

## Step 8 — Update the per-step call in `R/biting_process.R:204`

```r
kernels <- compute_atn_kernels(timestep, parameters, foim)
adult_mosquito_model_update(
  models[[s_i]]$.model,
  mu,
  foim,
  kernels$delta_atn,
  kernels$dn_atn,
  kernels$Lambda_i,
  kernels$rho_i,
  kernels$B_post,
  f
)
```

`foim` still routes only into `Sv[1]` (baseline); the ODE derivative handles
exposed rows via `Lambda_i[i]` for `i ≥ 2`.

---

## Step 9 — Update the EIR read-out (`R/biting_process.R:262`)

```r
calculate_infectious_compartmental <- function(solver_states, deltaq, spor_len) {
  iv_idx <- iv_block_indices(deltaq, spor_len)   # from Step 4
  max(sum(solver_states[iv_idx]), 0)
}
```

Pass `parameters$deltaq` and `parameters$spor_len` through
`calculate_infectious()` → `calculate_infectious_compartmental()`.

---

## Step 10 — Update compartmental rendering (`R/compartmental.R:95`)

`create_compartmental_rendering_process()` currently iterates over
`c(ODE_INDICES, ADULT_ODE_INDICES)`. Replace with:

```r
indices <- c(ODE_INDICES, make_adult_ode_indices(
  parameters$deltaq, parameters$spor_len))
```

Rendered column names will become `Sv1_species`, `Ev1_species`, … `Iv1_species`,
etc. Downstream scripts that depend on `Sm_`/`Pm_`/`Im_` column names will need
updating, but this is cosmetic and can be deferred post-verification.

---

## Step 11 — Run `Rcpp::compileAttributes()` + `devtools::load_all()`

After all C++ edits regenerate the Rcpp glue and confirm the package loads
cleanly with zero warnings.

---

## Step 12 — Verification (do not skip)

### 12a. Baseline (ATN off) regression test

With defaults (`Q0_atn = 0`, `n_atn = 1`, `deltaq = 1`, `spor_len = 1`):

- `delta_atn = 0` → the ODE collapses to `Sv[1]/Ev[1,1]/Iv[1]` only.
- `rho = 1/dem` and a single Erlang stage with `Ev_ratio = rho/(rho+mum)`.
- Compare `sum(Iv)`, total prevalence, and EIR against stock `malariasimulation`
  output. Expect match to ODE solver tolerance, not bit-identical (Erlang vs
  point delay differs by O(1/spor_len²) at `spor_len = 1`; increase `spor_len`
  to confirm convergence).

### 12b. ATN-on parity test

Match parameters (`deltaq`, `spor_len`, `gamma_atn`, `lambda_atn`, `p_atn`,
`t0_atn`, `Q0_atn`, kernels) between the fork and
`malariasimple_aitn_deterministic_v3.R` run via `dust2`. Compare `Svtot`,
`Evtot`, `Ivtot`, and EIR trajectories.

Add both as `tests/testthat/test-atn-mosquito.R`.

---

## Parameter / naming summary

| v3 symbol | malariasimulation parameter | Notes |
|---|---|---|
| `delayMos` | `parameters$dem` | EIP duration |
| `mu` | `parameters$mum[[i]]` | adult death rate |
| `Lambda` / `FOIvdel` | `foim` / `parameters$init_foim` | human→mosquito FOI |
| `mv0` | `m` (local in `initial_mosquito_counts`) | total adult density |
| `deltaqp1` | `parameters$deltaq + 1L` | computed, not stored |
| `rho` | `spor_len / dem` | computed in constructor and kernel function |
| `av` | bite rate `f` — passed into `adult_mosquito_model_update` | |
| `Svtot` | `sum(solver_states[sv_block_indices(...)])` | |
| `Ivtot` | `sum(solver_states[iv_block_indices(...)])` | replaces `Im` in EIR |
