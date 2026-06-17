# =====================================================================
# atn_convergence_check_par_v2.R — 2-D convergence (ATN only), PARALLEL
# =====================================================================
# Same as atn_convergence_check_v2.R but the ATN grid runs across cores.
#
# IMPORTANT (Windows / dev fork):
#  - Run from the repo ROOT, and run devtools::load_all() in THIS session
#    first (so the fork is compiled). Each worker is a fresh R process and
#    re-loads the dev fork with pkgload::load_all() (reusing that build) —
#    it does NOT pick up load_all() from the main session automatically.
#  - If worker loading ever fails (e.g. a DLL lock), either install the fork
#    once (devtools::install()) and swap the worker line to
#    library(malariasimulation), or fall back to the serial v2 script.
#
# Needs: ./dev/itn_params/{dat_res_cfp,dat_res_only}.rds
#        ./dev/atn_params/{eip_fit_hill_result,atn_fit_tra_main_bidir_splitnH}.rds
# =====================================================================

library(malariasimulation); library(dplyr); library(tidyr); library(ggplot2)
library(parallel)

fork_path <- getwd()   # repo root (run from here)

# ---- sweep settings ----
eir         <- 5
res         <- 0.4
atn_window  <- 10
deltaq_grid <- c(1L, 2L, 5L, 10L, 20L, 40L)   # exposed compartments
spor_grid   <- c(1L, 2L, 5L, 10L, 20L, 40L)   # infected (EIP) compartments

# ---- settings ----
n_steps                <- 13 * 365
old_distribution_times <- c(1, 4) * 365
new_distribution_times <- c(7, 10) * 365
retention_time         <- 588
old_covs <- rep(0.9, length(old_distribution_times))
new_covs <- rep(0.9, length(new_distribution_times))
distribution_times <- c(old_distribution_times, new_distribution_times)
n_old_dist <- length(old_distribution_times); n_new_dist <- length(new_distribution_times)
n_dist <- length(distribution_times)

# ---- median inputs (computed ONCE on main, then exported to workers) ----
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
rm(eip_samples, tra_samples)   # drop large Stan objects before exporting to workers

render_overrides <- list(
  human_population = 10000,
  prevalence_rendering_min_ages = 2*365,  prevalence_rendering_max_ages = 10*365,
  age_group_rendering_min_ages  = 2*365,  age_group_rendering_max_ages  = 10*365,
  clinical_incidence_rendering_min_ages = 0, clinical_incidence_rendering_max_ages = 100*365
)
form_overrides <- list(use_eip_hill = TRUE, use_bompard = TRUE)

# ---- run one arm at a given (deltaq, spor_len) ----
run_arm <- function(arm, eir, res, deltaq = 10L, atn_window = 10, spor_len = 10L) {
  cfp <- med_net(cfp_pars, res)
  struct <- list(deltaq = as.integer(deltaq), atn_window = atn_window,
                 spor_len = as.integer(spor_len))
  atn_overrides <- struct
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
    mutate(deltaq = deltaq, spor_len = spor_len,
           year = (timestep - new_distribution_times[1]) / 365,
           pfpr2to10 = n_detect_lm_730_3650 / n_age_730_3650,
           clin_inc  = n_inc_clinical_0_36500)
}

# ---- cluster ----
export_names <- ls()                                   # all small objects defined above
n_cores <- min(max(1L, detectCores() - 1L), length(deltaq_grid) * length(spor_grid))
message(sprintf("starting cluster on %d cores ...", n_cores))
cl <- makeCluster(n_cores)
clusterExport(cl, export_names)
invisible(clusterEvalQ(cl, {
  setwd(fork_path)
  suppressMessages(pkgload::load_all(fork_path, quiet = TRUE))  # load the DEV fork
  suppressMessages(library(dplyr))
}))
clusterSetRNGStream(cl, 123)   # reproducible parallel RNG (human IBM is stochastic)

# ---- run: baselines on main (fast), ATN grid in parallel ----
set.seed(123)
none_df <- lapply(spor_grid, function(sl)
  run_arm("none", eir, res, deltaq = 1L, atn_window = atn_window, spor_len = sl)
) |> bind_rows() |> filter(year >= 0, year <= 6)

combos     <- expand.grid(deltaq = deltaq_grid, spor_len = spor_grid)
combo_list <- split(combos, seq_len(nrow(combos)))
t0 <- proc.time()[["elapsed"]]
atn_list <- tryCatch(
  parLapply(cl, combo_list, function(co) {
    rt <- system.time(
      out <- run_arm("atn", eir, res, deltaq = co$deltaq, atn_window = atn_window, spor_len = co$spor_len)
    )[["elapsed"]]
    out$runtime_s <- rt
    out
  }),
  finally = stopCluster(cl)
)
atn_df <- bind_rows(atn_list) |> filter(year >= 0, year <= 6)
message(sprintf("grid complete in %.0f s on %d cores", proc.time()[["elapsed"]] - t0, n_cores))

# ---- plots (identical to serial v2) ----
atn_df  <- atn_df  |> mutate(deltaq_f = factor(deltaq, levels = deltaq_grid))
infect_lab <- labeller(spor_len = \(x) paste0("Infected compartments = ", x))

p_pfpr <- ggplot(atn_df, aes(year, pfpr2to10 * 100, colour = deltaq_f)) +
  geom_line(data = none_df, aes(year, pfpr2to10 * 100),
            colour = "grey70", linetype = "dashed", inherit.aes = FALSE) +
  geom_line(linewidth = 0.7) +
  facet_wrap(~ spor_len, labeller = infect_lab) +
  scale_colour_viridis_d(option = "C", end = 0.9) +
  labs(title = sprintf("ATN resolution convergence (atn_window = %g, EIR = %g, res = %g)",
                       atn_window, eir, res),
       subtitle = "coloured lines should converge as exposed compartments grow; grey = no future nets",
       x = "Years post future distribution",
       y = expression(italic(Pf) * PR[2-10] * " (%)"),
       colour = "Exposed\ncompartments") +
  theme_minimal(base_size = 13)

base_tot <- none_df |> group_by(spor_len) |>
  summarise(base_total = sum(clin_inc), .groups = "drop")
df_avert <- atn_df |> group_by(deltaq, deltaq_f, spor_len) |>
  summarise(arm_total = sum(clin_inc), .groups = "drop") |>
  left_join(base_tot, by = "spor_len") |>
  mutate(averted = base_total - arm_total)

p_avert <- ggplot(df_avert, aes(deltaq_f, averted, fill = deltaq_f)) +
  geom_col(width = 0.8) +
  facet_wrap(~ spor_len, labeller = infect_lab) +
  scale_fill_viridis_d(option = "C", end = 0.9) +
  labs(title = "Cases averted vs no future nets — should flatten as resolution grows",
       x = "Exposed compartments", y = "Clinical cases averted (summed over window)",
       fill = "Exposed\ncompartments") +
  theme_minimal(base_size = 13) + theme(legend.position = "none")

# ---- per-run compute time (s) over the deltaq x spor_len grid ----
# NB this is each combo's per-WORKER time; the wall-clock total is lower
# because runs ran concurrently.
runtimes <- atn_df |> distinct(deltaq, spor_len, runtime_s)
p_rt <- ggplot(runtimes, aes(factor(deltaq, levels = deltaq_grid),
                             factor(spor_len, levels = spor_grid), fill = runtime_s)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.1f", runtime_s)), size = 3) +
  scale_fill_viridis_c(option = "C") +
  coord_fixed() +
  labs(title = "Per-run compute time (s) — per worker",
       subtitle = "each combo's own runtime; wall-clock total is lower due to parallelism",
       x = "Exposed compartments (deltaq)",
       y = "Infected compartments (spor_len)", fill = "seconds") +
  theme_minimal(base_size = 13)

print(p_pfpr); print(p_avert); print(p_rt)
