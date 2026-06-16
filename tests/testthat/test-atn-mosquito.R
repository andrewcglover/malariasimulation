
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
    kernels <- compute_atn_kernels(t, parameters, parameters$init_foim)
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
    kernels <- compute_atn_kernels(t, parameters, parameters$init_foim)
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
