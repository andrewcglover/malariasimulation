# dev/generate-reference.R
#
# Generate regression reference data from the ORIGINAL (upstream) malariasimulation.
#
# When to run:
#   Run this script ONCE from a fresh R session in which the fork has NOT been
#   installed (devtools::load_all() does not count — callr spawns a clean
#   subprocess that only sees installed packages, so load_all() in the parent
#   session is fine).  If the fork has been installed over the original, first
#   reinstall the upstream:
#     remotes::install_github("mrc-ide/malariasimulation")
#   then run this script, then reinstall the fork.
#
# Output:
#   tests/testthat/fixtures/reference_sm_pm_im.rds  — list with mean EIR,
#   mean prevalence, and the parameter hash used, all from the original ODE.

if (!requireNamespace("callr", quietly = TRUE)) {
  stop("Install callr first:  install.packages('callr')")
}

out_path <- file.path("tests", "testthat", "fixtures", "reference_sm_pm_im.rds")
dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)

cat("Spawning clean subprocess to run original malariasimulation...\n")

ref <- callr::r(function() {
  library(malariasimulation)

  # --- confirm this is the original (no ATN params) ---
  p_check <- get_parameters()
  if (!is.null(p_check$deltaq)) {
    stop(
      "ATN fork appears to be installed. ",
      "Reinstall the original before generating reference data:\n",
      "  remotes::install_github('mrc-ide/malariasimulation')"
    )
  }

  set.seed(123L)
  p <- get_parameters(list(
    human_population                = 2000L,
    prevalence_rendering_min_ages   = 2L * 365L,
    prevalence_rendering_max_ages   = 10L * 365L
  ))
  p <- set_equilibrium(p, init_EIR = 10)
  sim <- run_simulation(365L * 3L, parameters = p)

  # Summarise over year 3 (post burn-in)
  yr3 <- sim[sim$timestep > 365L * 2L, ]
  prev <- yr3$n_detect_lm_730_3650 / yr3$n_730_3650

  list(
    mean_eir_yr3  = mean(yr3$EIR_gamb, na.rm = TRUE),
    mean_prev_yr3 = mean(prev, na.rm = TRUE),
    n_timesteps   = nrow(sim),
    pkg_version   = as.character(packageVersion("malariasimulation"))
  )
})

saveRDS(ref, out_path)

cat(sprintf("Reference saved to  %s\n", out_path))
cat(sprintf("  package version : %s\n",  ref$pkg_version))
cat(sprintf("  mean EIR yr3    : %.5f\n", ref$mean_eir_yr3))
cat(sprintf("  mean prev yr3   : %.5f\n", ref$mean_prev_yr3))
