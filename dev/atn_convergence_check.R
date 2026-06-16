# =====================================================================
# atn_convergence_check.R — exposure-axis resolution convergence test
# =====================================================================
# Checks that varying deltaq (number of time-since-exposure compartments)
# at a FIXED atn_window leaves the ATN results essentially unchanged —
# i.e. the generalisation samples the same continuous day-window more
# finely rather than changing the biology.
#
# Run AFTER devtools::load_all(), with working directory = repo root.
# Only the ATN arms are resolution-dependent; `none` is run once as the
# averted-cases baseline.
#
# Needs: ./dev/itn_params/{dat_res_cfp,dat_res_only}.rds
#        ./dev/atn_params/{eip_fit_hill_result,atn_fit_tra_main_bidir_splitnH}.rds
# =====================================================================

library(malariasimulation)
library(dplyr); library(tidyr); library(ggplot2)

# ---- what to sweep ----
eir          <- 5          # single transmission setting (convergence is setting-independent)
res          <- 0.4
atn_window   <- 10         # held fixed
deltaq_grid  <- c(10L, 20L, 40L)   # resolution: 10/20/40 (40 ~ quarter-day over 10 days)
spor_len_fix <- 10L

# ---- 1. settings (as in atn_local_test) ----
n_steps                <- 13 * 365
old_distribution_times <- c(1, 4) * 365
new_distribution_times <- c(7, 10) * 365
retention_time         <- 588
old_covs <- rep(0.9, length(old_distribution_times))
new_covs <- rep(0.9, length(new_distribution_times))
distribution_times <- c(old_distribution_times, new_distribution_times)
n_old_dist <- length(old_distribution_times); n_new_dist <- length(new_distribution_times)
n_dist <- length(distribution_times)

# ---- 2. median inputs ----
cfp_pars  <- readRDS("./dev/itn_params/dat_res_cfp.rds")  |> mutate(gamman = gamman * 365)
only_pars <- readRDS("./dev/itn_params/dat_res_only.rds") |> mutate(gamman = gamman * 365)
med_net <- function(pars, res) {
  out <- pars |> filter(resistance == res) |>
    summarise(dn0 = median(dn0), rn0 = median(rn0), gamman = median(gamman))
  if (nrow(out) == 0 || any(is.na(out))) stop("No net params for resistance ", res)
  out
}
eip_samples <- rstan::extract(readRDS("./dev/atn_params/eip_fit_hill_result.rds"))
tra_samples <- rstan::extract(readRDS("./dev/atn_params/atn_fit_tra_main_bidir_splitnH.rds"))
s_half_eip <- median(eip_samples$s_half); nH_eip <- median(eip_samples$nH)
rho_frac   <- median(eip_samples$r0_frac)
s_half_pre <- median(tra_samples$s_half_pre); nH_pre <- median(tra_samples$nH_pre)
B_max_post <- median(tra_samples$b_max)
s_half_post<- median(tra_samples$s_half_post); nH_post <- median(tra_samples$nH_post)
atn_halflife <- median(c(only_pars$gamman, cfp_pars$gamman))
gamma_atn    <- log(2) / atn_halflife

render_overrides <- list(
  human_population = 10000,
  prevalence_rendering_min_ages = 2*365,  prevalence_rendering_max_ages = 10*365,
  age_group_rendering_min_ages  = 2*365,  age_group_rendering_max_ages  = 10*365,
  clinical_incidence_rendering_min_ages = 0, clinical_incidence_rendering_max_ages = 100*365
)
form_overrides <- list(use_eip_hill = TRUE, use_bompard = TRUE)

# ---- 3. run one arm at a given resolution ----
# arm = "none" | "atn" | "pyr_atn"   (cfp omitted — not needed here)
run_arm <- function(arm, eir, res, deltaq = 10L, atn_window = 10, spor_len = 10L) {
  cfp <- med_net(cfp_pars, res); only <- med_net(only_pars, res)
  atn_overrides <- list()
  struct <- list(deltaq = as.integer(deltaq), atn_window = atn_window,
                 spor_len = as.integer(spor_len))

  if (arm == "none") {
    fut_cov <- rep(0, n_new_dist)
    fut_dn0 <- rep(cfp$dn0, n_new_dist); fut_rn <- rep(cfp$rn0, n_new_dist)
    fut_gam <- rep(cfp$gamman, n_new_dist); fut_rnm <- rep(0.24, n_new_dist)
  } else if (arm == "atn") {
    fut_cov <- new_covs
    fut_dn0 <- rep(0, n_new_dist); fut_rn <- rep(0.24, n_new_dist)
    fut_gam <- rep(atn_halflife, n_new_dist); fut_rnm <- rep(0.24 - 1e-9, n_new_dist)
    atn_overrides <- c(struct, list(
      p_atn = 0.9, gamma_atn = gamma_atn,
      s_half_eip = s_half_eip, nH_eip = nH_eip, rho_frac = rho_frac,
      s_half_pre = s_half_pre, nH_pre = nH_pre,
      B_max_post = B_max_post, s_half_post = s_half_post, nH_post = nH_post,
      Q0_atn = new_covs, t0_atn = new_distribution_times, n_atn = n_new_dist))
  } else if (arm == "pyr_atn") {
    fut_cov <- new_covs
    fut_dn0 <- rep(only$dn0, n_new_dist); fut_rn <- rep(only$rn0, n_new_dist)
    fut_gam <- rep(only$gamman, n_new_dist); fut_rnm <- rep(0.24, n_new_dist)
    atn_overrides <- c(struct, list(
      p_atn = 0.9, gamma_atn = gamma_atn,
      s_half_eip = s_half_eip, nH_eip = nH_eip, rho_frac = rho_frac,
      s_half_pre = s_half_pre, nH_pre = nH_pre,
      B_max_post = B_max_post, s_half_post = s_half_post, nH_post = nH_post,
      Q0_atn = new_covs, t0_atn = new_distribution_times, n_atn = n_new_dist,
      dn0_atn = only$dn0))
  }

  coverages <- c(old_covs, fut_cov)
  dn0_mat <- matrix(c(rep(cfp$dn0, n_old_dist), fut_dn0), ncol = 1)
  rn_mat  <- matrix(c(rep(cfp$rn0, n_old_dist), fut_rn),  ncol = 1)
  rnm_mat <- matrix(c(rep(0.24,    n_old_dist), fut_rnm), ncol = 1)
  gamman_v<- c(rep(cfp$gamman, n_old_dist), fut_gam)

  params <- get_parameters(overrides = c(render_overrides, form_overrides, atn_overrides)) |>
    set_bednets(timesteps = distribution_times, coverages = coverages,
                retention = retention_time,
                dn0 = dn0_mat, rn = rn_mat, rnm = rnm_mat, gamman = gamman_v) |>
    set_equilibrium(init_EIR = eir)
  r <- run_simulation(timesteps = n_steps, parameters = params)
  r |> as.data.frame() |>
    mutate(arm = arm, deltaq = deltaq,
           year = (timestep - new_distribution_times[1]) / 365,
           pfpr2to10 = n_detect_lm_730_3650 / n_age_730_3650,
           clin_inc  = n_inc_clinical_0_36500)
}

# ---- 4. run: none once + ATN arms across deltaq ----
message("baseline (none) ..."); none_df <- run_arm("none", eir, res) |>
  filter(year >= 0, year <= 6)

atn_df <- lapply(deltaq_grid, function(dq) {
  message(sprintf("ATN arms, deltaq = %d ...", dq))
  bind_rows(
    run_arm("atn",     eir, res, deltaq = dq, atn_window = atn_window, spor_len = spor_len_fix),
    run_arm("pyr_atn", eir, res, deltaq = dq, atn_window = atn_window, spor_len = spor_len_fix)
  )
}) |> bind_rows() |> filter(year >= 0, year <= 6)

arm_lab <- c(atn = "Future ATN", pyr_atn = "Future Pyr-ATN")
atn_df  <- atn_df  |> mutate(arm_f = factor(arm_lab[arm], levels = arm_lab))

# ---- 5. prevalence over time (lines = deltaq; grey dashed = no future nets) ----
p_pfpr <- ggplot(atn_df, aes(year, pfpr2to10 * 100, colour = factor(deltaq))) +
  geom_line(data = none_df,
            aes(year, pfpr2to10 * 100), colour = "grey60", linetype = "dashed",
            inherit.aes = FALSE) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~ arm_f) +
  labs(title = sprintf("Exposure-axis convergence (atn_window = %g, EIR = %g, res = %g)",
                       atn_window, eir, res),
       subtitle = "lines should overlap / tighten as deltaq grows; grey = no future nets",
       x = "Years post future distribution",
       y = expression(italic(Pf) * PR[2-10] * " (%)"), colour = "deltaq") +
  theme_minimal(base_size = 13)

# ---- 6. clinical cases averted vs no future nets ----
base_total <- sum(none_df$clin_inc)
df_avert <- atn_df |> group_by(arm_f, deltaq) |>
  summarise(arm_total = sum(clin_inc), .groups = "drop") |>
  mutate(averted = base_total - arm_total)

p_avert <- ggplot(df_avert, aes(factor(deltaq), averted, fill = factor(deltaq))) +
  geom_col(width = 0.7, alpha = 0.85) +
  facet_wrap(~ arm_f) +
  labs(title = "Cases averted vs no future nets — should be ~flat across deltaq",
       x = "deltaq", y = "Clinical cases averted (summed over window)", fill = "deltaq") +
  theme_minimal(base_size = 13) + theme(legend.position = "none")

print(p_pfpr); print(p_avert)

# ---- 7. independence: deltaq and spor_len set together (no dimension error) ----
chk <- run_arm("atn", eir, res, deltaq = 20L, atn_window = atn_window, spor_len = 15L)
message(sprintf("independence check (deltaq=20, spor_len=15): pfpr range %.3f–%.3f",
                min(chk$pfpr2to10), max(chk$pfpr2to10)))
