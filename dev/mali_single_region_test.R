# dev/mali_single_region_test.R
#
# Smoke test: run all 4 arms sequentially for the highest-EIR Mali admin-1
# region, then print a summary table.
#
# Usage (from the fork root, after devtools::load_all()):
#   source("dev/mali_single_region_test.R")

# Source the run script for its functions + data; skip the cluster launch.
options(mali_test_mode = TRUE)
source("dev/mali_projection_run.R")
options(mali_test_mode = NULL)

# ── Pick highest-EIR region ───────────────────────────────────────────────────
eir_by_region <- vapply(regions, function(rg) {
  site_row <- site_obj$sites[site_obj$sites$name_1 == rg, , drop = FALSE]
  ms       <- site::subset_site(site_obj, site_row)
  ms$eir$eir
}, numeric(1))

test_region <- names(which.max(eir_by_region))
message(sprintf("\nHighest-EIR region: %s  (EIR = %.1f)\n", test_region, max(eir_by_region)))

# ── Run all 4 arms for that region ───────────────────────────────────────────
results <- lapply(arms, function(arm) {
  message(sprintf("  arm = %-8s ...", arm), appendLF = FALSE)
  t0  <- proc.time()[["elapsed"]]
  out <- run_one(list(region = test_region, arm = arm))
  dt  <- proc.time()[["elapsed"]] - t0
  message(sprintf(" done (%.0f s)", dt))
  out
})

df <- dplyr::bind_rows(results)

# ── Summary: mean future PfPR and clin incidence per arm ─────────────────────
future_df <- df[df$year_rel >= 0, ]
summary_tbl <- future_df |>
  dplyr::group_by(arm) |>
  dplyr::summarise(
    mean_pfpr2to10  = round(mean(pfpr2to10,  na.rm = TRUE), 4),
    mean_clin_inc   = round(mean(clin_inc,   na.rm = TRUE), 1),
    .groups = "drop"
  )

message("\nFuture-period summary (years 0–6 relative to ", future_yr0, "):")
print(as.data.frame(summary_tbl))
