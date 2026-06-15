# Regression test: fork with ATN off vs original Sm/Pm/Im malariasimulation.
#
# The fork replaces the 3-compartment Sm/Pm/Im ODE (fixed EIP delay via deque)
# with a 2-D Sv/Ev/Iv Erlang chain (spor_len stages).  With ATN off and
# spor_len = 10 (default), Erlang EIP approximates the fixed delay; long-run
# EIR and prevalence should agree with the original to within ~10%.
#
# Fixture: tests/testthat/fixtures/reference_sm_pm_im.rds
#   Generate once with the original package:  source("dev/generate-reference.R")

test_that('fork (ATN off) EIR matches original Sm/Pm/Im to within 10%', {
  fixture <- testthat::test_path("fixtures", "reference_sm_pm_im.rds")
  skip_if_not(
    file.exists(fixture),
    "Reference fixture missing — run dev/generate-reference.R with the original package"
  )
  ref <- readRDS(fixture)

  set.seed(123L)
  p <- get_parameters(list(
    human_population              = 2000L,
    prevalence_rendering_min_ages = 2L * 365L,
    prevalence_rendering_max_ages = 10L * 365L
  ))
  p <- set_equilibrium(p, init_EIR = 10)
  sim <- run_simulation(365L * 3L, parameters = p)

  yr3     <- sim[sim$timestep > 365L * 2L, ]
  mean_eir <- mean(yr3$EIR_gamb, na.rm = TRUE)

  expect_equal(mean_eir, ref$mean_eir_yr3, tolerance = 0.10,
    label = "mean EIR (year 3)")
})

test_that('fork (ATN off) prevalence matches original Sm/Pm/Im to within 10%', {
  fixture <- testthat::test_path("fixtures", "reference_sm_pm_im.rds")
  skip_if_not(
    file.exists(fixture),
    "Reference fixture missing — run dev/generate-reference.R with the original package"
  )
  ref <- readRDS(fixture)

  set.seed(123L)
  p <- get_parameters(list(
    human_population              = 2000L,
    prevalence_rendering_min_ages = 2L * 365L,
    prevalence_rendering_max_ages = 10L * 365L
  ))
  p <- set_equilibrium(p, init_EIR = 10)
  sim <- run_simulation(365L * 3L, parameters = p)

  yr3  <- sim[sim$timestep > 365L * 2L, ]
  prev <- yr3$n_detect_lm_730_3650 / yr3$n_730_3650
  mean_prev <- mean(prev, na.rm = TRUE)

  expect_equal(mean_prev, ref$mean_prev_yr3, tolerance = 0.10,
    label = "mean prevalence 2-10yr (year 3)")
})
