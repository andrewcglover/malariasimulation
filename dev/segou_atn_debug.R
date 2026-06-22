# =====================================================================
# segou_atn_debug.R   (DIAGNOSTIC — parallel PSOCK cluster)
#
#   Goal: localize why Ségou shows the lowest ATN-exposed mosquito fraction
#   of any Mali region in the pure `atn` arm. Renders the exact per-species
#   inputs to the ATN exposure pathway (dbg_*) so each region's exposed
#   fraction can be reconciled against av_da = a * delta_atn * contact_factor.
#
#   Run from the FORK ROOT:
#       source("dev/segou_atn_debug.R")
# =====================================================================

setwd("C:/Users/ag4218/Local/GitHub/malariasimulation")
options(country_test_mode = TRUE)
source("dev/c24med_projection_run.R")   # defines build_params + all globals

DEBUG_POP <- 20000L
render_overrides$human_population <- DEBUG_POP

debug_regions <- c("Ségou", "Sikasso", "Gao", "Mopti")
debug_arms    <- c("atn", "pyr_cfp")

grid_df <- expand.grid(region = debug_regions, arm = debug_arms,
                       stringsAsFactors = FALSE)
rows    <- split(grid_df, seq_len(nrow(grid_df)))
n_sims  <- length(rows)

# Cores: min(#sims, 16, detectCores()-2). Generalisable for larger future runs.
n_cores <- min(n_sims, 16L, max(1L, parallel::detectCores() - 2L))
message(sprintf("Debug sweep: %d sims (%d regions x %d arms) on %d cores, pop=%d",
                n_sims, length(debug_regions), length(debug_arms), n_cores, DEBUG_POP))

debug_run_one <- function(row) {
  region <- row$region; arm <- row$arm
  t0 <- proc.time()[["elapsed"]]
  tryCatch({
    p <- build_params(region, arm)
    p$atn_debug <- TRUE
    r <- run_simulation(timesteps = n_steps, parameters = p)
    d <- as.data.frame(r)
    d$region   <- region
    d$arm      <- arm
    d$year_rel <- (d$timestep - future_start_day) / 365
    message(sprintf("[%s x %s] done in %.0f s", region, arm,
                    proc.time()[["elapsed"]] - t0))
    d
  }, error = function(e) {
    list(.__error__ = TRUE, region = region, arm = arm, message = conditionMessage(e))
  })
}

# setup_strategy="sequential": workers connect one at a time (avoids parallel-setup
# race condition that causes "invalid connection" in Rscript sessions on Windows).
# No outfile arg: default discards worker stdout (outfile="" pipe-conflicts with
# the Rscript stdout redirect and kills connections after workers start).
cl <- parallel::makeCluster(n_cores, setup_strategy = "sequential")
on.exit(parallel::stopCluster(cl), add = TRUE)

parallel::clusterExport(cl, "FORK_PATH")
parallel::clusterEvalQ(cl, {
  suppressMessages({ library(site); library(dplyr); library(tidyr) })
  pkgload::load_all(FORK_PATH, quiet = TRUE, compile = FALSE)
})
parallel::clusterExport(cl, c(
  "build_params", "build_past_schedule", "build_future_schedule", "med_net",
  "expand_interventions", "debug_run_one",
  "site_obj", "regions", "start_year", "hist_last", "future_yr0",
  "future_start_day", "future_campaign_days", "n_steps", "n_future_years",
  "cfp_pars", "only_pars", "pbo_pars", "atn_kern", "gamma_atn", "ANTIMAL_HL_DAYS",
  "cd_cov", "campaign_cov", "cd_interval", "retention_time", "deltaq_use",
  "render_overrides", "form_overrides", "human_pop",
  "future_interventions"
))

t0_total <- proc.time()[["elapsed"]]
res_list <- parallel::parLapplyLB(cl, rows, debug_run_one)
message(sprintf("All %d sims done in %.0f s", n_sims,
                proc.time()[["elapsed"]] - t0_total))

is_error <- vapply(res_list, function(x) isTRUE(x$.__error__), logical(1))
if (any(is_error)) {
  message(sprintf("WARNING: %d job(s) failed:", sum(is_error)))
  print(dplyr::bind_rows(lapply(res_list[is_error], function(x)
    data.frame(region = x$region, arm = x$arm, message = x$message))))
}
dbg <- dplyr::bind_rows(res_list[!is_error])

out <- "dev/outputs/mli_c24med_segou_atn_debug.rds"
saveRDS(list(
  results = dbg,
  meta    = list(human_pop = DEBUG_POP, future_yr0 = future_yr0,
                 n_future_years = n_future_years,
                 regions = debug_regions, arms = debug_arms)
), out)
message("Saved: ", out)

# ---------------------------------------------------------------------
# Summary (future window): is av_da region-uniform within each arm?
# ---------------------------------------------------------------------
fy <- n_future_years
for (sp in c("gambiae", "arabiensis", "funestus")) {
  s <- dbg |>
    dplyr::filter(year_rel >= 0, year_rel <= fy) |>
    dplyr::group_by(arm, region) |>
    dplyr::summarise(
      a           = mean(.data[[paste0("dbg_a_", sp)]]),
      W           = mean(.data[[paste0("dbg_W_", sp)]]),
      Z           = mean(.data[[paste0("dbg_Z_", sp)]]),
      delta_atn   = mean(.data[[paste0("dbg_delta_atn_", sp)]]),
      contact_fac = mean(.data[[paste0("dbg_contact_factor_", sp)]]),
      av_da       = mean(.data[[paste0("dbg_av_da_", sp)]]),
      Sv_exp_frac = sum(.data[[paste0("Sv_exposed_", sp, "_count")]]) /
                    sum(.data[[paste0("Sv_unexposed_", sp, "_count")]] +
                        .data[[paste0("Sv_exposed_", sp, "_count")]]),
      .groups = "drop"
    ) |>
    dplyr::arrange(arm, region)
  cat("\n=== ", sp, " (future-window means) ===\n", sep = "")
  print(as.data.frame(s), digits = 4)
}
