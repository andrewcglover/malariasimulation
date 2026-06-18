
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
