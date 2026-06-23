
# ── Step 12 verification tests for the ATN adult-mosquito ODE ─────────────────

# ── 12a: Baseline (ATN off) — ODE equilibrium is stable ──────────────────────

test_that('Adult ODE stays at Erlang equilibrium with ATN off', {
  parameters <- get_parameters()
  parameters <- set_equilibrium(parameters, 50.)
  f  <- parameters$blood_meal_rates[[1]]
  mu <- parameters$mum[[1]]

  models  <- parameterise_mosquito_models(parameters, 365L)
  solvers <- parameterise_solvers(models, parameters)

  equilibrium <- initial_mosquito_counts(
    parameters, 1L, parameters$init_foim, parameters$total_M
  )

  # Adult-only indices: aquatic (E/L/P) drift ~2e-4 from size_t truncation of
  # total_M in the aquatic ODE sub-model — original malariasimulation behaviour.
  adult_idx <- c(
    sv_block_indices(parameters$deltaq),
    ev_block_indices(parameters$deltaq, parameters$spor_len),
    iv_block_indices(parameters$deltaq, parameters$spor_len)
  )

  # Initial adult state must match Erlang equilibrium exactly
  expect_equal(solvers[[1]]$get_states()[adult_idx], equilibrium[adult_idx],
               tolerance = 1e-8)

  # Run for one year; adult compartments must stay at equilibrium
  for (t in seq_len(365L)) {
    kernels <- compute_atn_kernels(t, parameters, parameters$init_foim, 1L)
    adult_mosquito_model_update(
      models[[1]]$.model,
      mu,
      parameters$init_foim,
      f * kernels$delta_atn,
      kernels$delta_atn,
      kernels$dn_atn,
      kernels$Lambda0_t,
      kernels$Lambda_i,
      kernels$rho_i,
      kernels$B_post,
      f
    )
    solvers[[1]]$step()
  }

  expect_equal(solvers[[1]]$get_states()[adult_idx], equilibrium[adult_idx],
               tolerance = 1e-4)
})


test_that('calculate_infectious_compartmental returns positive Ivtot at equilibrium', {
  parameters <- get_parameters()
  parameters <- set_equilibrium(parameters, 50.)
  eq <- initial_mosquito_counts(
    parameters, 1L, parameters$init_foim, parameters$total_M
  )
  ivtot <- calculate_infectious_compartmental(eq, parameters)
  expect_gt(ivtot, 0)
})


# ── 12b: ATN on — extra mortality reduces infectious mosquitoes ───────────────

test_that('ATN extra mortality reduces total infectious mosquitoes', {
  # Use small compartment counts for speed; effect must be visible after 1 year
  parameters <- get_parameters(list(
    deltaq   = 3L,
    spor_len = 5L,
    # ATN active from day 1 with high extra mortality
    p_atn    = 0.8,
    Q0_atn   = 0.9,
    t0_atn   = 1L,
    lambda_atn = 0,      # no coverage decay
    gamma_atn  = 0,      # no drug-effect decay
    dn0_atn    = 0.5     # 50% extra daily mortality for exposed mosquitoes
  ))
  parameters <- set_equilibrium(parameters, 50.)
  f  <- parameters$blood_meal_rates[[1]]
  mu <- parameters$mum[[1]]

  models  <- parameterise_mosquito_models(parameters, 365L)
  solvers <- parameterise_solvers(models, parameters)

  iv_idx <- iv_block_indices(parameters$deltaq, parameters$spor_len)
  initial_Ivtot <- sum(solvers[[1]]$get_states()[iv_idx])

  for (t in seq_len(365L)) {
    kernels <- compute_atn_kernels(t, parameters, parameters$init_foim, 1L)
    adult_mosquito_model_update(
      models[[1]]$.model,
      mu,
      parameters$init_foim,
      f * kernels$delta_atn,
      kernels$delta_atn,
      kernels$dn_atn,
      kernels$Lambda0_t,
      kernels$Lambda_i,
      kernels$rho_i,
      kernels$B_post,
      f
    )
    solvers[[1]]$step()
  }

  final_Ivtot <- sum(solvers[[1]]$get_states()[iv_idx])
  # ATN kills mosquitoes — infectious compartment should have dropped
  expect_lt(final_Ivtot, initial_Ivtot * 0.9)
})


# ── 12d: compute_atn_kernels invariants — pure R, no ODE stepping ─────────────
#
# Each test calls compute_atn_kernels() directly; no solver is stepped so these
# run in milliseconds and are safe to select-all-and-run after devtools::load_all().

test_that('compute_atn_kernels: ATN-off defaults collapse to baseline', {
  # With p_atn=0 and Q0_atn=0 (defaults), coverage Q_t=0 and delta_atn=0.
  # Lambda_i should equal foim everywhere; rho_i should equal rho everywhere.
  parameters <- get_parameters()
  parameters <- set_equilibrium(parameters, 50.)
  foim <- parameters$init_foim
  rho  <- parameters$spor_len / parameters$dem
  deltaqp1 <- parameters$deltaq + 1L

  k <- compute_atn_kernels(1L, parameters, foim, 1L)

  expect_equal(k$delta_atn, 0)
  expect_equal(k$dn_atn,    0)
  expect_equal(k$Lambda0_t, foim)
  expect_equal(k$Lambda_i,  rep(foim, deltaqp1))
  expect_equal(k$rho_i,     rep(rho,  deltaqp1))
})

test_that('compute_atn_kernels: zero exposure before first distribution round', {
  # Q_t=0 when timestep < t0_atn, regardless of p_atn / Q0_atn values.
  parameters <- get_parameters(list(
    p_atn      = 0.9,
    Q0_atn     = 0.8,
    t0_atn     = 200L,
    lambda_atn = 0       # no decay — makes Q_t unambiguous after t0_atn
  ))
  parameters <- set_equilibrium(parameters, 50.)
  foim <- parameters$init_foim

  k <- compute_atn_kernels(199L, parameters, foim, 1L)

  expect_equal(k$delta_atn, 0)
  expect_equal(k$dn_atn,    0)
  expect_equal(k$Lambda0_t, foim)
})

test_that('compute_atn_kernels: p_atn=0 gates mosquito exposure even with coverage', {
  # Q_t > 0 after t0_atn, but p_atn=0 means no mosquito makes drug contact.
  parameters <- get_parameters(list(
    p_atn      = 0,
    Q0_atn     = 0.8,
    t0_atn     = 1L,
    lambda_atn = 0
  ))
  parameters <- set_equilibrium(parameters, 50.)
  foim <- parameters$init_foim

  k <- compute_atn_kernels(100L, parameters, foim, 1L)

  expect_equal(k$delta_atn, 0)
})

test_that('compute_atn_kernels: active kernel gives correct delta_atn', {
  # With no coverage decay (lambda_atn=0) and a single round of Q0_atn=0.8,
  # Q_t == 0.8 at all timesteps >= t0_atn, so delta_atn = p_atn * phi * 0.8.
  parameters <- get_parameters(list(
    p_atn      = 0.9,
    Q0_atn     = 0.8,
    t0_atn     = 1L,
    lambda_atn = 0
  ))
  parameters <- set_equilibrium(parameters, 50.)
  foim <- parameters$init_foim
  phi  <- parameters$phi_bednets[[1]]

  k <- compute_atn_kernels(100L, parameters, foim, 1L)

  expect_equal(k$delta_atn, 0.9 * phi * 0.8, tolerance = 1e-10)
  expect_gt(k$delta_atn, 0)
})

test_that('compute_atn_kernels: minimum dimensions (deltaq=1, spor_len=1) are finite', {
  # Guards the s[1]<-0 NaN-prevention for the Hill-kernel computation
  # and the spor_len=1 / deltaq=1 empty-loop edge cases.
  parameters <- get_parameters(list(
    deltaq     = 1L,
    spor_len   = 1L,
    p_atn      = 0.9,
    Q0_atn     = 0.8,
    t0_atn     = 1L,
    lambda_atn = 0
  ))
  parameters <- set_equilibrium(parameters, 50.)
  foim <- parameters$init_foim

  k <- compute_atn_kernels(50L, parameters, foim, 1L)

  expect_equal(length(k$Lambda_i), 2L)   # deltaq + 1 = 2
  expect_equal(length(k$rho_i),    2L)
  expect_equal(length(k$B_post),   1L)   # spor_len = 1
  expect_true(all(is.finite(unlist(k))))
})

test_that('compute_atn_kernels: delta_atn is species-specific via phi_bednets', {
  # Concern 3: phi_bednets[[species]] must index correctly for each species.
  # With three species, each should give p_atn * phi_s * Q_t for its own phi.
  parameters <- get_parameters(list(
    p_atn      = 0.9,
    Q0_atn     = 0.8,
    t0_atn     = 1L,
    lambda_atn = 0
  ))
  parameters <- set_species(
    parameters,
    list(gamb_params, arab_params, fun_params),
    c(0.4, 0.3, 0.3)
  )
  parameters <- set_equilibrium(parameters, 50.)
  foim <- parameters$init_foim

  for (s in 1:3) {
    k   <- compute_atn_kernels(100L, parameters, foim, s)
    phi <- parameters$phi_bednets[[s]]
    expect_equal(k$delta_atn, 0.9 * phi * 0.8, tolerance = 1e-10,
                 label = paste0('species ', s))
  }
})

# ── 12c: Full simulation with ATN parameters completes without error ──────────

test_that('run_simulation completes with ATN parameters set', {
  parameters <- get_parameters(list(
    human_population = 500L,
    deltaq   = 2L,
    spor_len = 5L,
    p_atn    = 0.5,
    Q0_atn   = 0.8,
    t0_atn   = 100L,
    lambda_atn = 0.003,
    gamma_atn  = 0.003,
    B_max_post = 0.3,    # peak pre- & post-blocking = 0.3; Lambda00 = foim*(1-0.3) derived
    dn0_atn    = 0.2,
    prevalence_rendering_min_ages = 2 * 365,
    prevalence_rendering_max_ages = 10 * 365
  ))
  parameters <- set_equilibrium(parameters, 20.)
  sim <- run_simulation(200L, parameters = parameters)
  expect_equal(nrow(sim), 200L)
  expect_true('EIR_gamb' %in% names(sim))
})

# ── 12e: lambda_atn auto-derive from set_bednets retention ───────────────────

test_that('set_bednets auto-derives lambda_atn = 1/retention when lambda_atn is NULL', {
  parameters <- get_parameters()
  expect_null(parameters$lambda_atn)   # default is NULL (auto-derive sentinel)
  n_sp <- length(parameters$species)
  parameters <- set_bednets(
    parameters,
    timesteps = 100L,
    coverages = 0.5,
    retention = 2000,
    dn0 = matrix(0,   nrow = 1, ncol = n_sp),
    rn  = matrix(0.24, nrow = 1, ncol = n_sp),
    rnm = matrix(0.24 - 1e-9, nrow = 1, ncol = n_sp),
    gamman = 365
  )
  expect_equal(parameters$lambda_atn, 1 / 2000, tolerance = 1e-12)
})

test_that('set_bednets respects an explicitly set lambda_atn (does not override)', {
  parameters <- get_parameters(list(lambda_atn = 0))   # explicit 0
  n_sp <- length(parameters$species)
  parameters <- set_bednets(
    parameters,
    timesteps = 100L,
    coverages = 0.5,
    retention = 2000,
    dn0 = matrix(0,   nrow = 1, ncol = n_sp),
    rn  = matrix(0.24, nrow = 1, ncol = n_sp),
    rnm = matrix(0.24 - 1e-9, nrow = 1, ncol = n_sp),
    gamman = 365
  )
  expect_equal(parameters$lambda_atn, 0)   # must NOT be overwritten
})

test_that('compute_atn_kernels falls back to lambda_atn=0 when no set_bednets called', {
  # NULL lambda_atn (no set_bednets) must not error and must give no coverage decay
  parameters <- get_parameters(list(
    p_atn  = 0.9,
    Q0_atn = 0.8,
    t0_atn = 1L
    # lambda_atn left as NULL
  ))
  parameters <- set_equilibrium(parameters, 50.)
  foim <- parameters$init_foim
  # If lambda_atn=0 fallback is working, coverage at t=365 equals Q0_atn (no decay)
  k <- compute_atn_kernels(365L, parameters, foim, 1L)
  expect_equal(k$delta_atn, 0.9 * parameters$phi_bednets[[1]] * 0.8, tolerance = 1e-10)
})

test_that('set_bednets with logistic retention warns and sets lambda_atn = 1/half_life', {
  parameters <- get_parameters()
  n_sp <- length(parameters$species)
  expect_warning(
    parameters <- set_bednets(
      parameters,
      timesteps = 100L,
      coverages = 0.5,
      logistic_half_life = 1500,
      logistic_k = 20,
      dn0 = matrix(0,   nrow = 1, ncol = n_sp),
      rn  = matrix(0.24, nrow = 1, ncol = n_sp),
      rnm = matrix(0.24 - 1e-9, nrow = 1, ncol = n_sp),
      gamman = 365
    ),
    regexp = "logistic net retention"
  )
  expect_equal(parameters$lambda_atn, 1 / 1500, tolerance = 1e-12)
})

# ── 12f: contact_factor — barrier-repelled mosquitoes get dosed ───────────────
#
# contact_factor = (sn + rnm) / (1 - rnm) scales av_da to include mosquitoes that
# physically touch the net (barrier-repelled, prob rnm) but are not fed-and-survived.
# Only chemical excito-repellency (rn - rnm) prevents net contact entirely.
# Denominator (1 - rnm) is the untreated-net floor: fixed ~0.76, never collapses.
# Numerator (sn + rnm) excludes pyrethroid-killed mosquitoes (dn); they cannot transmit.
# All tests call compute_atn_kernels() directly with set_bednets set up so that
# t0_atn matches a bednet schedule row (as in the Mali pipeline).

test_that('compute_atn_kernels: contact_factor = 1 when no set_bednets called', {
  # NULL parameters$bednet_timesteps -> fallback to 1 (no adjustment).
  parameters <- get_parameters(list(
    p_atn = 0.9, Q0_atn = 0.8, t0_atn = 1L, lambda_atn = 0
  ))
  parameters <- set_equilibrium(parameters, 50.)
  k <- compute_atn_kernels(100L, parameters, parameters$init_foim, 1L)
  expect_equal(k$contact_factor, 1)
})

test_that('compute_atn_kernels: contact_factor = 1/(1-rnm) for non-insecticidal ATN', {
  # rn0 = rnm, dn0 = 0: sn = 1-rnm (constant); (sn+rnm)/sn = 1/(1-rnm).
  rnm_val <- 0.24 - 1e-9
  parameters <- get_parameters(list(
    p_atn = 0.9, Q0_atn = 0.8, t0_atn = 100L, lambda_atn = 0
  ))
  n_sp <- length(parameters$species)
  parameters <- set_bednets(
    parameters,
    timesteps = 100L, coverages = 0.8, retention = 5000,
    dn0    = matrix(0,       nrow = 1, ncol = n_sp),
    rn     = matrix(0.24,    nrow = 1, ncol = n_sp),
    rnm    = matrix(rnm_val, nrow = 1, ncol = n_sp),
    gamman = 365 * 5
  )
  parameters <- set_equilibrium(parameters, 50.)
  k <- compute_atn_kernels(200L, parameters, parameters$init_foim, 1L)
  expect_equal(k$contact_factor, 1 / (1 - rnm_val), tolerance = 1e-9)
})

test_that('compute_atn_kernels: contact_factor for fresh Pyr-ATN matches (sn+rnm)/(1-rnm) at dt=0', {
  # At dt=0, rn(0)=rn0, dn(0)=dn0, sn(0)=1-rn0-dn0.
  # contact_factor = (sn+rnm)/(1-rnm): bounded, excludes killed (dn), preserves coupling.
  # Must be < 1/(1-rnm) (the non-insecticidal ATN value): Pyr-ATN < ATN ordering.
  rn0_val <- 0.5; rnm_val <- 0.24; dn0_val <- 0.3
  parameters <- get_parameters(list(
    p_atn = 0.9, Q0_atn = 0.8, t0_atn = 100L, lambda_atn = 0
  ))
  n_sp <- length(parameters$species)
  parameters <- set_bednets(
    parameters,
    timesteps = 100L, coverages = 0.8, retention = 5000,
    dn0    = matrix(dn0_val, nrow = 1, ncol = n_sp),
    rn     = matrix(rn0_val, nrow = 1, ncol = n_sp),
    rnm    = matrix(rnm_val, nrow = 1, ncol = n_sp),
    gamman = 365 * 5   # slow decay -> dt=0 approximation exact at timestep=t0
  )
  parameters <- set_equilibrium(parameters, 50.)
  k <- compute_atn_kernels(100L, parameters, parameters$init_foim, 1L)
  sn_expected <- 1 - rn0_val - dn0_val
  cf_expected <- (sn_expected + rnm_val) / (1 - rnm_val)   # bounded denominator
  expect_equal(k$contact_factor, cf_expected, tolerance = 1e-9)
  # Pyr-ATN contact_factor < non-insecticidal ATN value (ordering check)
  expect_lt(k$contact_factor, 1 / (1 - rnm_val))
  # factor < 1: pyrethroid removes more contacts than barrier-repelled adds back
  expect_lt(k$contact_factor, 1)
})

test_that('compute_atn_kernels: contact_factor decays toward 1/(1-rnm) as Pyr-ATN ages', {
  # With small gamman, rn->rnm and dn->0 quickly, so contact_factor -> 1/(1-rnm).
  rnm_val <- 0.24
  parameters <- get_parameters(list(
    p_atn = 0.9, Q0_atn = 0.8, t0_atn = 100L, lambda_atn = 0
  ))
  n_sp <- length(parameters$species)
  parameters <- set_bednets(
    parameters,
    timesteps = 100L, coverages = 0.8, retention = 5000,
    dn0    = matrix(0.3,     nrow = 1, ncol = n_sp),
    rn     = matrix(0.5,     nrow = 1, ncol = n_sp),
    rnm    = matrix(rnm_val, nrow = 1, ncol = n_sp),
    gamman = 10   # fast decay: at dt = 500 >> gamman, rn~rnm, dn~0
  )
  parameters <- set_equilibrium(parameters, 50.)
  k <- compute_atn_kernels(600L, parameters, parameters$init_foim, 1L)
  expect_equal(k$contact_factor, 1 / (1 - rnm_val), tolerance = 1e-6)
})

# ── §12g: displacement events (atn_displace_t0/Q0) ───────────────────────────

test_that('compute_atn_kernels: empty displacement vectors leave existing arms unchanged', {
  # Backward-compatibility: atn_displace_t0/Q0 default to numeric(0).
  # A pyr_cfp_atn-style config with no displacement should give exactly the
  # same repl_factor as before the parameter was introduced.
  parameters <- get_parameters(list(
    p_atn = 0.9, n_atn = 2L,
    t0_atn  = c(100L, 400L),
    Q0_atn  = c(0.8, 0.3),
    lambda_atn = 0
    # atn_displace_t0/Q0 intentionally omitted -> default numeric(0)
  ))
  parameters <- set_equilibrium(parameters, 50.)
  # Before either drug event fires: Q_t = 0, delta_atn = 0
  k_before <- compute_atn_kernels(50L, parameters, parameters$init_foim, 1L)
  expect_equal(k_before$delta_atn, 0)
  # Between the two drug events: only event 1 active; no displacement yet
  k_mid <- compute_atn_kernels(300L, parameters, parameters$init_foim, 1L)
  expect_equal(k_mid$delta_atn,
               parameters$p_atn * parameters$phi_bednets[[1]] * 0.8,
               tolerance = 1e-9)
  # After both drug events: Q_t = Q0[1]*(1-Q0[2]) + Q0[2] (no displacement)
  k_after <- compute_atn_kernels(500L, parameters, parameters$init_foim, 1L)
  Q_expected <- 0.8 * (1 - 0.3) + 0.3
  expect_equal(k_after$delta_atn,
               parameters$p_atn * parameters$phi_bednets[[1]] * Q_expected,
               tolerance = 1e-9)
})

test_that('compute_atn_kernels: displacement event collapses Q_atn_t correctly', {
  # One ATN CD event (t=100, Q0=0.8) followed by one non-drug displacement (t=300, Q0=0.9).
  # Before t=300:  repl_factor for event 1 = 1 (no later events fired yet).
  # After  t=300:  repl_factor for event 1 = (1 - 0.9) = 0.1.
  parameters <- get_parameters(list(
    p_atn           = 0.9,
    n_atn           = 1L,
    t0_atn          = 100L,
    Q0_atn          = 0.8,
    lambda_atn      = 0,
    atn_displace_t0 = 300L,
    atn_displace_Q0 = 0.9
  ))
  parameters <- set_equilibrium(parameters, 50.)
  phi <- parameters$phi_bednets[[1]]

  # Timestep 200: drug event active, displacement not yet fired -> Q_t = 0.8
  k_pre <- compute_atn_kernels(200L, parameters, parameters$init_foim, 1L)
  expect_equal(k_pre$delta_atn, parameters$p_atn * phi * 0.8, tolerance = 1e-9)

  # Timestep 400: displacement has fired -> repl_factor[1] = (1 - 0.9) = 0.1
  k_post <- compute_atn_kernels(400L, parameters, parameters$init_foim, 1L)
  expect_equal(k_post$delta_atn, parameters$p_atn * phi * 0.8 * 0.1, tolerance = 1e-9)
})

test_that('compute_atn_kernels: CD event after displacement is unaffected by that displacement', {
  # Displacement at t=200 must NOT reduce a drug CD event that fires at t=300
  # (displacement is chronologically earlier than the drug event).
  parameters <- get_parameters(list(
    p_atn           = 0.9,
    n_atn           = 1L,
    t0_atn          = 300L,
    Q0_atn          = 0.5,
    lambda_atn      = 0,
    atn_displace_t0 = 200L,
    atn_displace_Q0 = 0.9
  ))
  parameters <- set_equilibrium(parameters, 50.)
  phi <- parameters$phi_bednets[[1]]

  # Timestep 400: drug event active; displacement at t=200 is EARLIER, so not in
  # "later" set for event 1 -> repl_factor[1] = 1 -> Q_t = 0.5
  k <- compute_atn_kernels(400L, parameters, parameters$init_foim, 1L)
  expect_equal(k$delta_atn, parameters$p_atn * phi * 0.5, tolerance = 1e-9)
})
