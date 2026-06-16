# =====================================================================
# atn_local_test.R  —  stripped-down local sanity test of the ATN fork
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
n_steps                <- 13 * 365          # 13 years, daily timesteps
old_distribution_times <- c(1, 4) * 365     # past nets (Pyr-CFP in every arm)
new_distribution_times <- c(7, 10) * 365    # future campaign rounds (2 rounds)
retention_time         <- 588

old_coverage <- 0.9
new_coverage <- 0.9

eir_levels <- 5#c(5, 50)
res_levels <- 0.4#c(0.4, 0.9)

# ATN structural resolution. deltaq = exposure compartments (fork default is 1;
# malariasimple used 10). With atn_window at its default of 10, deltaq = 10 -> kappa = 1.
# spor_len (EIP stages) already defaults to 10, so it isn't set here.
deltaq_use <- 10L

# Derived
distribution_times <- c(old_distribution_times, new_distribution_times)
n_old_dist <- length(old_distribution_times)
n_new_dist <- length(new_distribution_times)
n_dist     <- length(distribution_times)
old_covs   <- rep(old_coverage, n_old_dist)
new_covs   <- rep(new_coverage, n_new_dist)
all_covs   <- c(old_covs, new_covs)

# ---------------------------------------------------------------------
# 2. Median input parameters
# ---------------------------------------------------------------------
# Net efficacy (Pyr-CFP, Pyr-only) by resistance — median across draws.
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

# ATN kernels — median of Stan posteriors.
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

# ATN drug-effect half-life -> gamma_atn. v6 reuses net-retention half-lives
# (gamman, in days after *365) as the ATN half-life distribution; take the median.
# NB confirm the gamman -> half-life semantics/units match your intent.
atn_halflife <- median(c(only_pars$gamman, cfp_pars$gamman))
gamma_atn    <- log(2) / atn_halflife

# ---------------------------------------------------------------------
# 3. Output rendering config
# ---------------------------------------------------------------------
render_overrides <- list(
  human_population                      = 10000,     # smoother single-run incidence
  prevalence_rendering_min_ages         = 2 * 365,
  prevalence_rendering_max_ages         = 10 * 365,
  age_group_rendering_min_ages          = 2 * 365,   # denominator (n_age_730_3650)
  age_group_rendering_max_ages          = 10 * 365,
  clinical_incidence_rendering_min_ages = 0,
  clinical_incidence_rendering_max_ages = 100 * 365
)

# Hill + Bompard: both TRUE (kept explicit as a harmless safety net).
form_overrides <- list(use_eip_hill = TRUE, use_bompard = TRUE)

# ---------------------------------------------------------------------
# 4. Run one arm of one (EIR, resistance) cell
#    arm = "none" | "cfp" | "atn" | "pyr_atn"
# ---------------------------------------------------------------------
run_arm <- function(arm, eir, res) {

  cfp  <- med_net(cfp_pars,  res)
  only <- med_net(only_pars, res)

  # future-net efficacy + coverage + ATN overrides by arm (past = Pyr-CFP always)
  atn_overrides <- list()
  if (arm == "none") {
    fut_cov <- rep(0, n_new_dist)
    fut_dn0 <- rep(cfp$dn0, n_new_dist); fut_rn <- rep(cfp$rn0, n_new_dist)
    fut_gam <- rep(cfp$gamman, n_new_dist)
    fut_rnm <- rep(0.24, n_new_dist)
  } else if (arm == "cfp") {
    fut_cov <- new_covs
    fut_dn0 <- rep(cfp$dn0, n_new_dist); fut_rn <- rep(cfp$rn0, n_new_dist)
    fut_gam <- rep(cfp$gamman, n_new_dist)
    fut_rnm <- rep(0.24, n_new_dist)
  } else if (arm == "atn") {                       # antimalarial net, no insecticide
    fut_cov <- new_covs
    fut_dn0 <- rep(0, n_new_dist); fut_rn <- rep(0.24, n_new_dist)
    fut_gam <- rep(atn_halflife, n_new_dist)
    fut_rnm <- rep(0.24 - 1e-9, n_new_dist)  # no insecticide -> rnm ≈ rn; strict < required
    atn_overrides <- list(
      p_atn = 0.9, deltaq = deltaq_use,
      gamma_atn = gamma_atn,
      s_half_eip = s_half_eip, nH_eip = nH_eip, rho_frac = rho_frac,
      s_half_pre = s_half_pre, nH_pre = nH_pre,
      B_max_post = B_max_post, s_half_post = s_half_post, nH_post = nH_post,
      Q0_atn = new_covs, t0_atn = new_distribution_times, n_atn = n_new_dist
    )
  } else if (arm == "pyr_atn") {                   # pyrethroid net + antimalarial
    fut_cov <- new_covs
    fut_dn0 <- rep(only$dn0, n_new_dist); fut_rn <- rep(only$rn0, n_new_dist)
    fut_gam <- rep(only$gamman, n_new_dist)
    fut_rnm <- rep(0.24, n_new_dist)
    atn_overrides <- list(
      p_atn = 0.9, deltaq = deltaq_use,
      gamma_atn = gamma_atn,
      s_half_eip = s_half_eip, nH_eip = nH_eip, rho_frac = rho_frac,
      s_half_pre = s_half_pre, nH_pre = nH_pre,
      B_max_post = B_max_post, s_half_post = s_half_post, nH_post = nH_post,
      Q0_atn = new_covs, t0_atn = new_distribution_times, n_atn = n_new_dist,
      dn0_atn = only$dn0
    )
  }

  coverages <- c(old_covs, fut_cov)
  dn0_mat   <- matrix(c(rep(cfp$dn0, n_old_dist), fut_dn0), ncol = 1)  # rows=dists, cols=species
  rn_mat    <- matrix(c(rep(cfp$rn0, n_old_dist), fut_rn),  ncol = 1)
  rnm_mat   <- matrix(c(rep(0.24, n_old_dist), fut_rnm), ncol = 1)
  gamman_v  <- c(rep(cfp$gamman, n_old_dist), fut_gam)

  params <- get_parameters(overrides = c(render_overrides, form_overrides, atn_overrides))

  params <- params |>
    set_bednets(
      timesteps = distribution_times,
      coverages = coverages,
      retention = retention_time,
      dn0 = dn0_mat, rn = rn_mat, rnm = rnm_mat, gamman = gamman_v
    ) |>
    set_equilibrium(init_EIR = eir)

  r <- run_simulation(timesteps = n_steps, parameters = params)

  r |> as.data.frame() |>
    mutate(
      arm = arm, EIR = eir, resistance = res,
      year      = (timestep - new_distribution_times[1]) / 365,
      pfpr2to10 = n_detect_lm_730_3650 / n_age_730_3650,   # microscopy PfPR2-10 (v6 used n_detect)
      clin_inc  = n_inc_clinical_0_36500
    )
}

# ---------------------------------------------------------------------
# 5. Run the grid x 4 arms
# ---------------------------------------------------------------------
arms <- c("none", "cfp", "atn", "pyr_atn")
grid <- expand.grid(eir = eir_levels, res = res_levels, arm = arms,
                    stringsAsFactors = FALSE)

n_runs   <- nrow(grid)
t_start  <- proc.time()[["elapsed"]]

run_arm_progress <- function(arm, eir, res, i) {
  message(sprintf("[%2d/%d]  arm=%-7s  EIR=%2d  res=%2.0f%%  ...",
                  i, n_runs, arm, eir, res * 100))
  t0     <- proc.time()[["elapsed"]]
  result <- run_arm(arm, eir, res)
  elapsed <- proc.time()[["elapsed"]] - t0
  elapsed_total <- proc.time()[["elapsed"]] - t_start
  message(sprintf("        done in %.1f s  (total %.1f s)", elapsed, elapsed_total))
  result
}

df_full <- Map(run_arm_progress, grid$arm, grid$eir, grid$res, seq_len(n_runs)) |>
  bind_rows() |>
  filter(year >= 0, year <= 6)

message(sprintf("All %d simulations complete in %.1f s.", n_runs,
                proc.time()[["elapsed"]] - t_start))

# ---------------------------------------------------------------------
# 6. Plotting helpers
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

print(p_pfpr)
print(p_avert)
