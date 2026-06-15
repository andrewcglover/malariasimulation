test_that('compartmental always gives positive infectious', {
  parameters <- get_parameters()
  state_len <- 3L + (parameters$deltaq + 1L) * (2L + parameters$spor_len)
  solver_states <- rep(0, state_len)
  iv_idx <- iv_block_indices(parameters$deltaq, parameters$spor_len)
  solver_states[[iv_idx[[1]]]] <- -1e-10
  expect_gte(calculate_infectious_compartmental(solver_states, parameters), 0)
})

test_that('gonotrophic_cycle cannot be negative', {
  params <- get_parameters()
  vparams <- gamb_params
  vparams$blood_meal_rates <- 5
  expect_error(set_species(params, list(vparams), 1))
})
