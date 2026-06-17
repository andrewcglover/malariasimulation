# =====================================================================
# atn_local_test_cd_v1.R
#   Local sanity test of the ATN fork, EXTENDED with continuous
#   distribution (CD) of nets alongside the 3-yearly mass campaigns.
#   Synthetic single site (get_parameters + set_bednets + set_equilibrium),
#   four arms: none / cfp / atn / pyr_atn.
#
#   Lineage: atn_local_test_v2.R. Only the net SCHEDULE changes:
#   each campaign is now backed by monthly mini top-up distributions that
#   hold a maintained "CD floor" of net use between campaigns.
# =====================================================================
# Run AFTER devtools::load_all(), with working directory = repo root.
#
# Needs input files present:
#   ./dev/itn_params/dat_res_cfp.rds, dat_res_only.rds            (ITN net efficacy)
#   ./dev/atn_params/eip_fit_hill_result.rds                      (Stan: EIP fit)
#   ./dev/atn_params/atn_fit_tra_main_bidir_splitnH.rds           (Stan: TRA fit)
# =====================================================================

library(malariasimulation)
library(dplyr)
library(tidyr)
library(ggplot2)

# ---------------------------------------------------------------------
# 1. Settings  (median inputs, no draws)
# ---------------------------------------------------------------------
n_steps        <- 13 * 365          # 13 years, daily timesteps
retention_time <- 588               # mean net-RETENTION time (days). Drives net-USE decay.

# Mass-campaign schedule: every 36 months. Past (Pyr-CFP, common to all arms)
# at years 1 & 4; future (arm-specific) at years 7 & 10.
old_campaign_times <- c(1, 4)  * 365
new_campaign_times <- c(7, 10) * 365
future_start       <- new_campaign_times[1]      # arms diverge from here

campaign_cov <- 0.9                 # campaign COVERAGE: fraction who receive a new net
                                    # at a mass campaign (random allocation). NB this is
                                    # the input to set_bednets; realised post-campaign
                                    # USE (U0, below) is a touch higher.

# --- Continuous distribution (CD) ------------------------------------
# A fraction CD_frac of post-campaign net USE originates from CD; that sets a
# "floor" of use that monthly top-ups maintain between campaigns.
CD_frac     <- 0.155                # ~15.5% (95% CI 13.5-17.4) from the prior paper
cd_interval <- 365 / 12             # monthly top-up cadence (~30.4 days)

eir_levels <- 5     # c(5, 50)
res_levels <- 0.4   # c(0.4, 0.9)

# ATN structural resolution (see atn_local_test_v2.R).
deltaq_use <- 10L

# ---------------------------------------------------------------------
# 1b. CD coverage math  (random / proportional net replacement)
# ---------------------------------------------------------------------
# Under random/proportional allocation a distribution of coverage c overlaps the
# current users at random, so the UNCOVERED fraction multiplies:
#     1 - use_after = (1 - c) * (1 - use_before)   <=>   use_after = c + (1-c)*use_before
#
# (a) FLOOR. The CD floor is 15.5% of the post-campaign USE U0 (not of the 0.9
#     coverage). U0 itself = c + (1-c)*u_pre, so it depends on the pre-campaign use,
#     which depends on the floor -> mildly circular. Taking u_pre = the maintained
#     floor (programme-intrinsic, independent of campaign cadence) closes the loop:
#         floor = CD_frac*(c + (1-c)*floor)
#       => floor = CD_frac*c / (1 - CD_frac*(1-c))
#     (For c=0.9: floor=0.1417, U0=0.9142, CD_frac*U0=floor. Checked.)
#
# (b) MONTHLY TOP-UP. A single locked-in coverage `cd_cov` whose decay-then-topup
#     steady state equals the floor:
#         floor = cd_cov / (1 - (1-cd_cov)*d)   =>   cd_cov = floor*(1-d)/(1-floor*d)
#     (For the above floor: cd_cov=0.00825; steady state reproduces 0.1417.)
d_cd     <- exp(-cd_interval / retention_time)
cd_floor <- CD_frac * campaign_cov / (1 - CD_frac * (1 - campaign_cov))
U0       <- campaign_cov + (1 - campaign_cov) * cd_floor   # implied post-campaign use
cd_cov   <- cd_floor * (1 - d_cd) / (1 - cd_floor * d_cd)

message(sprintf(
  "Campaign coverage = %.3f -> post-campaign use U0 = %.4f | CD floor = %.4f (= %.1f%% of U0) | monthly top-up coverage = %.4f (%.2f%%)",
  campaign_cov, U0, cd_floor, 100 * CD_frac, cd_cov, 100 * cd_cov))

# ---------------------------------------------------------------------
# 2. Median input parameters
# ---------------------------------------------------------------------
message("Loading ITN net-efficacy parameters...")
cfp_pars  <- readRDS("./dev/itn_params/dat_res_cfp.rds")  |> mutate(gamman = gamman * 365)
only_pars <- readRDS("./dev/itn_params/dat_res_only.rds") |> mutate(gamman = gamman * 365)

med_net <- function(pars, res) {
  out <- pars |> filter(resistance == res) |>
    summarise(dn0 = median(dn0), rn0 = median(rn0), gamman = median(gamman))
  if (nrow(out) == 0 || any(is.na(out)))
    stop(sprintf("No net params for resistance = %s — check the RDS resistance grid.", res))
  out
}

message("Loading ATN Stan posteriors...")
eip_samples <- rstan::extract(readRDS("./dev/atn_params/eip_fit_hill_result.rds"))
tra_samples <- rstan::extract(readRDS("./dev/atn_params/atn_fit_tra_main_bidir_splitnH.rds"))

s_half_eip  <- median(eip_samples$s_half)
nH_eip      <- median(eip_samples$nH)
rho_frac    <- median(eip_samples$r0_frac)
s_half_pre  <- median(tra_samples$s_half_pre)
nH_pre      <- median(tra_samples$nH_pre)
B_max_post  <- median(tra_samples$b_max)
s_half_post <- median(tra_samples$s_half_post)
nH_post     <- median(tra_samples$nH_post)

# ATN drug-effect half-life -> gamma_atn (reuse net-retention half-lives, *365).
atn_halflife <- median(c(only_pars$gamman, cfp_pars$gamman))
gamma_atn    <- log(2) / atn_halflife

# ---------------------------------------------------------------------
# 3. Output rendering config
# ---------------------------------------------------------------------
render_overrides <- list(
  human_population                      = 10000,
  prevalence_rendering_min_ages         = 2 * 365,
  prevalence_rendering_max_ages         = 10 * 365,
  age_group_rendering_min_ages          = 2 * 365,
  age_group_rendering_max_ages          = 10 * 365,
  clinical_incidence_rendering_min_ages = 0,
  clinical_incidence_rendering_max_ages = 100 * 365
)

form_overrides <- list(use_eip_hill = TRUE, use_bompard = TRUE)

# ---------------------------------------------------------------------
# 4. Build the CD net schedule for one arm
# ---------------------------------------------------------------------
# Returns a list of per-distribution vectors ready for set_bednets(), plus the
# ATN event vectors (future ATN top-ups only). The monthly grid lands EXACTLY on
# campaign times (36 months = 36 * 365/12 = 1095 days), so campaigns are just the
# grid points that coincide with a campaign date.
#
# Net type by era:
#   - t < future_start                : Pyr-CFP (common past), every arm.
#   - t >= future_start               : arm-specific future net type.
#   - arm == "none"                   : NO future distributions at all.
build_schedule <- function(arm, cfp, only) {

  # Monthly grid from the first campaign to the end of the run (integer days).
  grid <- unique(round(seq(old_campaign_times[1], n_steps, by = cd_interval)))
  campaign_times <- c(old_campaign_times, new_campaign_times)

  is_camp <- vapply(grid, function(t) any(abs(t - campaign_times) < 1), logical(1))
  is_past <- grid <  future_start
  is_fut  <- grid >= future_start

  # "none": drop every future distribution (no future nets).
  if (arm == "none") {
    keep <- is_past
    grid <- grid[keep]; is_camp <- is_camp[keep]
    is_past <- rep(TRUE, length(grid)); is_fut <- rep(FALSE, length(grid))
  }

  n <- length(grid)
  coverages <- ifelse(is_camp, campaign_cov, cd_cov)

  # Per-distribution net-efficacy params. Default = Pyr-CFP (used for all past
  # distributions and for the cfp future arm).
  dn0    <- rep(cfp$dn0,    n)
  rn     <- rep(cfp$rn0,    n)
  rnm    <- rep(0.24,       n)
  gamman <- rep(cfp$gamman, n)

  # Overlay the arm-specific FUTURE net type.
  if (arm == "atn") {                 # antimalarial net, no insecticide
    dn0[is_fut]    <- 0
    rn[is_fut]     <- 0.24
    rnm[is_fut]    <- 0.24 - 1e-9     # no insecticide -> rnm ~ rn; strict < required
    gamman[is_fut] <- atn_halflife
  } else if (arm == "pyr_atn") {      # pyrethroid + antimalarial
    dn0[is_fut]    <- only$dn0
    rn[is_fut]     <- only$rn0
    rnm[is_fut]    <- 0.24
    gamman[is_fut] <- only$gamman
  }
  # arm %in% c("none","cfp") -> future (if any) stays Pyr-CFP defaults.

  # ATN event vectors: only the FUTURE distributions deliver antimalarial.
  # Each future top-up is an ATN event; Q0_atn carries its distribution coverage.
  fut_times <- grid[is_fut]
  fut_cov   <- coverages[is_fut]

  list(
    timesteps = grid,
    coverages = coverages,
    dn0_mat   = matrix(dn0,    ncol = 1),   # rows = distributions, cols = species
    rn_mat    = matrix(rn,     ncol = 1),
    rnm_mat   = matrix(rnm,    ncol = 1),
    gamman_v  = gamman,
    atn_times = fut_times,
    atn_cov   = fut_cov,
    n_atn     = length(fut_times)
  )
}

# ---------------------------------------------------------------------
# 5. Run one arm of one (EIR, resistance) cell
#    arm = "none" | "cfp" | "atn" | "pyr_atn"
# ---------------------------------------------------------------------
run_arm <- function(arm, eir, res) {

  cfp  <- med_net(cfp_pars,  res)
  only <- med_net(only_pars, res)

  sch <- build_schedule(arm, cfp, only)

  # ATN overrides (only the antimalarial arms switch the mechanism on).
  atn_overrides <- list()
  if (arm == "atn") {
    atn_overrides <- list(
      p_atn = 0.9, deltaq = deltaq_use, gamma_atn = gamma_atn,
      s_half_eip = s_half_eip, nH_eip = nH_eip, rho_frac = rho_frac,
      s_half_pre = s_half_pre, nH_pre = nH_pre,
      B_max_post = B_max_post, s_half_post = s_half_post, nH_post = nH_post,
      Q0_atn = sch$atn_cov, t0_atn = sch$atn_times, n_atn = sch$n_atn
    )
  } else if (arm == "pyr_atn") {
    atn_overrides <- list(
      p_atn = 0.9, deltaq = deltaq_use, gamma_atn = gamma_atn,
      s_half_eip = s_half_eip, nH_eip = nH_eip, rho_frac = rho_frac,
      s_half_pre = s_half_pre, nH_pre = nH_pre,
      B_max_post = B_max_post, s_half_post = s_half_post, nH_post = nH_post,
      Q0_atn = sch$atn_cov, t0_atn = sch$atn_times, n_atn = sch$n_atn,
      dn0_atn = only$dn0
    )
  }

  params <- get_parameters(overrides = c(render_overrides, form_overrides, atn_overrides))

  params <- params |>
    set_bednets(
      timesteps = sch$timesteps,
      coverages = sch$coverages,
      retention = retention_time,
      dn0 = sch$dn0_mat, rn = sch$rn_mat, rnm = sch$rnm_mat, gamman = sch$gamman_v
    ) |>
    set_equilibrium(init_EIR = eir)

  r <- run_simulation(timesteps = n_steps, parameters = params)

  r |> as.data.frame() |>
    mutate(
      arm = arm, EIR = eir, resistance = res,
      year      = (timestep - future_start) / 365,
      pfpr2to10 = n_detect_lm_730_3650 / n_age_730_3650,
      clin_inc  = n_inc_clinical_0_36500
    )
}

# ---------------------------------------------------------------------
# 6. Run the grid x 4 arms
# ---------------------------------------------------------------------
arms <- c("none", "cfp", "atn", "pyr_atn")
grid_cells <- expand.grid(eir = eir_levels, res = res_levels, arm = arms,
                          stringsAsFactors = FALSE)

n_runs  <- nrow(grid_cells)
t_start <- proc.time()[["elapsed"]]

run_arm_progress <- function(arm, eir, res, i) {
  message(sprintf("[%2d/%d]  arm=%-7s  EIR=%2d  res=%2.0f%%  ...",
                  i, n_runs, arm, eir, res * 100))
  t0     <- proc.time()[["elapsed"]]
  result <- run_arm(arm, eir, res)
  message(sprintf("        done in %.1f s  (total %.1f s)",
                  proc.time()[["elapsed"]] - t0, proc.time()[["elapsed"]] - t_start))
  result
}

df_full <- Map(run_arm_progress, grid_cells$arm, grid_cells$eir, grid_cells$res,
               seq_len(n_runs)) |>
  bind_rows() |>
  filter(year >= 0, year <= 6)

message(sprintf("All %d simulations complete in %.1f s.", n_runs,
                proc.time()[["elapsed"]] - t_start))

# ---------------------------------------------------------------------
# 7. Plotting helpers
# ---------------------------------------------------------------------
message("Building plots...")
arm_labels <- c(none = "No future nets", cfp = "Future Pyr-CFP",
                atn = "Future ATN", pyr_atn = "Future Pyr-ATN")
arm_cols   <- c("No future nets" = "grey50", "Future Pyr-CFP" = "#009988",
                "Future ATN" = "#EE7733", "Future Pyr-ATN" = "#CC3311")
facet_lab  <- labeller(
  EIR        = \(x) paste0("EIR = ", x),
  resistance = \(x) paste0("Resistance = ", as.numeric(x) * 100, "%")
)

df_full <- df_full |> mutate(arm_f = factor(arm_labels[arm], levels = arm_labels))

# Prevalence over time
p_pfpr <- ggplot(df_full, aes(year, pfpr2to10 * 100, colour = arm_f)) +
  geom_line(linewidth = 0.7) +
  scale_colour_manual(values = arm_cols) +
  facet_grid(EIR ~ resistance, labeller = facet_lab, scales = "free_y") +
  theme_minimal(base_size = 13) +
  labs(x = "Years post future distribution",
       y = expression(italic(Pf) * PR[2-10] * " (%)"), colour = "")

# Clinical cases averted vs "no future nets" (summed over the window)
base_inc <- df_full |> filter(arm == "none") |>
  group_by(EIR, resistance) |> summarise(base_total = sum(clin_inc), .groups = "drop")

df_avert <- df_full |> filter(arm != "none") |>
  group_by(arm, EIR, resistance) |>
  summarise(arm_total = sum(clin_inc), .groups = "drop") |>
  left_join(base_inc, by = c("EIR", "resistance")) |>
  mutate(averted = base_total - arm_total,
         arm_f = factor(arm_labels[arm], levels = arm_labels))

p_avert <- ggplot(df_avert, aes(arm_f, averted, fill = arm_f)) +
  geom_col(width = 0.7, alpha = 0.85) +
  scale_fill_manual(values = arm_cols) +
  facet_grid(EIR ~ resistance, labeller = facet_lab, scales = "free_y") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
  labs(x = "", y = "Clinical cases averted vs no future nets (summed over window)")

# ---------------------------------------------------------------------
# 8. CD floor check — realised net use vs the proportional-replacement
#    prediction and the target floor.
# ---------------------------------------------------------------------
# `predict_use_daily()` reconstructs the EXPECTED daily net-use trajectory
# implied by the schedule under the proportional-replacement assumption
# (decay one day, then use_after = c + (1-c)*use on distribution days). If
# set_bednets shares that semantics, realised use should track this closely;
# a large gap flags a different replacement rule.
predict_use_daily <- function(sch, t_max) {
  cov_at <- setNames(sch$coverages, as.character(sch$timesteps))
  u <- numeric(t_max); cur <- 0
  for (t in seq_len(t_max)) {
    cur <- cur * exp(-1 / retention_time)            # one day of net-use decay
    key <- as.character(t)
    if (!is.na(cov_at[key])) { c <- cov_at[[key]]; cur <- c + (1 - c) * cur }
    u[t] <- cur
  }
  data.frame(timestep = seq_len(t_max), pred_use = u)
}

usage_col <- intersect(c("n_use_net", "net_usage"), names(df_full))
if (length(usage_col) >= 1) {

  # Realised net use (auto-detect count vs fraction).
  raw_max <- max(df_full[[usage_col[1]]], na.rm = TRUE)
  to_frac <- if (raw_max > 1.5) render_overrides$human_population else 1
  df_use <- df_full |>
    mutate(net_use = .data[[usage_col[1]]] / to_frac, year_abs = timestep / 365)

  # Predicted trajectory per arm (schedule depends only on arm, not EIR/res).
  pred <- lapply(arms, function(a) {
    sch <- build_schedule(a, med_net(cfp_pars, res_levels[1]), med_net(only_pars, res_levels[1]))
    p <- predict_use_daily(sch, n_steps); p$arm <- a; p
  }) |> bind_rows()

  df_use <- df_use |> left_join(pred, by = c("arm", "timestep"))

  # Numeric check over the future window (years 0-6).
  cd_check <- df_use |>
    group_by(arm, EIR, resistance) |>
    summarise(
      real_mean   = mean(net_use),  pred_mean   = mean(pred_use),
      real_trough = min(net_use),   pred_trough = min(pred_use),
      mean_abs_diff = mean(abs(net_use - pred_use)),
      max_abs_diff  = max(abs(net_use - pred_use)),
      .groups = "drop"
    ) |>
    mutate(target_floor = cd_floor)

  message("\n--- CD floor / proportional-replacement check (future window, yrs 0-6) ---")
  print(as.data.frame(cd_check), row.names = FALSE, digits = 4)

  tol <- 0.02   # mean abs diff tolerance (allows ~1-day timing offset on the sawtooth)
  bad <- cd_check |> filter(arm != "none", mean_abs_diff > tol)
  if (nrow(bad) > 0) {
    warning(sprintf(
      "Realised net use diverges from the proportional-replacement prediction (mean abs diff > %.3f) in %d CD arm(s). set_bednets may use a different replacement rule — recheck the cd_cov mapping.",
      tol, nrow(bad)))
  } else {
    message(sprintf("OK: realised net use tracks the proportional-replacement prediction (mean abs diff <= %.3f) in all CD arms.", tol))
  }

  # Usage plot: realised (solid) vs predicted (dashed), with the floor line.
  p_use <- ggplot(df_use, aes(year_abs, net_use, colour = arm_f)) +
    geom_line(linewidth = 0.6) +
    geom_line(aes(y = pred_use), linetype = "dotted", linewidth = 0.5) +
    geom_hline(yintercept = cd_floor, linetype = "dashed", colour = "grey40") +
    scale_colour_manual(values = arm_cols) +
    facet_grid(EIR ~ resistance, labeller = facet_lab) +
    theme_minimal(base_size = 13) +
    labs(x = "Year (absolute)", y = "Net usage (fraction)", colour = "",
         caption = "Solid = realised, dotted = proportional-replacement prediction, dashed = CD floor")
  print(p_use)
} else {
  message("Note: no net-usage column found in output — skipping CD-floor check/plot.")
}

print(p_pfpr)
print(p_avert)
