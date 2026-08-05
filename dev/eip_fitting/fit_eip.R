# fit_eip.R
# Driver for eip_fit.stan — fits rho_0 and zeta (impaired-EIP recovery
# parameters) against paired sporozoite-positive binomial counts, using
# rstan and a shared OLRE per experimental row.

suppressPackageStartupMessages({
  library(rstan)
  library(posterior)
  library(bayesplot)
  library(dplyr)
  library(ggplot2)
  library(jsonlite)
})

# rstan options
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

set.seed(2025)

# ---- Paths ----
here_dir <- tryCatch({
  args <- commandArgs(trailingOnly = FALSE)
  fn <- sub("--file=", "", args[grep("--file=", args)])
  if (length(fn) && nzchar(fn)) {
    dirname(normalizePath(fn))
  } else if (requireNamespace("rstudioapi", quietly = TRUE) &&
             rstudioapi::isAvailable()) {
    dirname(rstudioapi::getActiveDocumentContext()$path)
  } else {
    getwd()
  }
}, error = function(e) getwd())
setwd(here_dir)

data_path  <- "./exp_data/SPZ_pooled_summary_APR26.csv"
model_path <- "eip_fit.stan"

# ---- Priors (single source of truth) ----
source("priors_atn_main.R")
prior_data <- prior_stan_data(priors)

# ---- Data ----
dat <- read.csv(data_path)

stopifnot(
  all(c("exp_id", "post_time", "day",
        "n_CTL", "pos_CTL", "n_200mg", "pos_200mg") %in% names(dat)),
  all(dat$pos_CTL   <= dat$n_CTL),
  all(dat$pos_200mg <= dat$n_200mg)
  # No sign constraint on post_time: negative values mean exposure
  # before infection and are handled by the Delta < 0 branch in
  # eip_fit.stan's likelihood.
)

# Re-index exp_id to 1..K in case gaps exist
dat$exp_idx <- as.integer(factor(dat$exp_id))

stan_data <- c(
  list(
    N       = nrow(dat),
    K       = length(unique(dat$exp_idx)),
    exp_id  = dat$exp_idx,
    delta   = dat$post_time / 24,
    tau_obs = dat$day,
    n_ctl   = dat$n_CTL,
    pos_ctl = dat$pos_CTL,
    n_trt   = dat$n_200mg,
    pos_trt = dat$pos_200mg,
    Delta_r = 10L
  ),
  prior_data
)

# ---- Compile ----
mod <- stan_model(file = model_path)

# ---- Sample ----
# Note: rstan's `iter` is TOTAL per chain including warmup.
# iter = 3000, warmup = 1000 -> 2000 sampling iters per chain.
fit <- sampling(
  mod,
  data    = stan_data,
  seed    = 2025,
  chains  = 4,
  iter    = 3000,
  warmup  = 1000,
  control = list(adapt_delta = 0.95, max_treedepth = 12),
  refresh = 500
)

# ---- Diagnostics ----
cat("\n---- HMC diagnostics ----\n")
check_hmc_diagnostics(fit)

cat("\n---- Key parameter summary ----\n")
key_vars <- c("log_eip_excess", "eip_excess", "eip_baseline",
              "rho", "rho0", "r0_frac",
              "zeta", "log_zeta", "zeta_halflife",
              "sigma_exp")
print(summary(fit, pars = key_vars)$summary)

cat("\n---- OLRE offsets (beta_exp) ----\n")
print(summary(fit, pars = "beta_exp")$summary)

# ---- Pairs plot for the ridge diagnostics ----
# Use the unconstrained primary parameters so the geometry HMC sees is visible.
post_arr <- as.array(fit, pars = c("log_eip_excess", "r0_frac", "log_zeta", "sigma_exp"))
p_pairs <- mcmc_pairs(post_arr,
                      off_diag_args = list(size = 0.4, alpha = 0.25))
print(p_pairs)
ggsave("eip_pairs.png", p_pairs, width = 9, height = 9, dpi = 200)

# ---- Trace plot ----
p_trace <- mcmc_trace(post_arr)
ggsave("eip_trace.png", p_trace, width = 10, height = 6, dpi = 200)

# ---- Posterior predictive checks ----
pp_ctl <- rstan::extract(fit, pars = "pos_ctl_rep")$pos_ctl_rep
pp_trt <- rstan::extract(fit, pars = "pos_trt_rep")$pos_trt_rep

ppc_ctl <- ppc_intervals(
  y    = dat$pos_CTL,
  yrep = pp_ctl,
  x    = seq_len(nrow(dat))
) +
  labs(title    = "Control arm: posterior predictive",
       subtitle = "x = row index (= exp_id after re-indexing)",
       x = "row", y = "sporozoite-positive count")
ggsave("ppc_ctl.png", ppc_ctl, width = 9, height = 5, dpi = 200)

ppc_trt <- ppc_intervals(
  y    = dat$pos_200mg,
  yrep = pp_trt,
  x    = seq_len(nrow(dat))
) +
  labs(title    = "Treatment (200 mg) arm: posterior predictive",
       subtitle = "x = row index (= exp_id after re-indexing)",
       x = "row", y = "sporozoite-positive count")
ggsave("ppc_trt.png", ppc_trt, width = 9, height = 5, dpi = 200)

# ---- Save fit ----
# NB: avoid "eip_fit.rds" — that filename collides with rstan's auto_write
# cache for eip_fit.stan, and the next stan_model() call would try to
# deserialise the stanfit as a compiled model and fail with a cryptic
# "no method for coercing this S4 class to a vector" error.
saveRDS(fit, file = "eip_fit_result.rds")
write_priors_json("eip_fit_priors.json", priors)
cat("\nSaved: eip_fit_result.rds, eip_fit_priors.json, eip_pairs.png, eip_trace.png, ppc_ctl.png, ppc_trt.png\n")

# ---- Identifiability check ----
post <- rstan::extract(fit, pars = c("rho", "rho0", "zeta", "r0_frac"))

cat("\n---- Joint posterior correlations ----\n")
cor_mat <- cor(cbind(rho0 = post$rho0, zeta = post$zeta, r0_frac = post$r0_frac))
print(cor_mat)

# Effective cumulative rate at readout for Delta = +3 d, tau_obs = 10 and 13.
# These are the identifiable combinations; report alongside individual params.
A_post <- (post$rho - post$rho0) / post$zeta
s10_3  <- post$rho * 10 - A_post * (1 - exp(-post$zeta * 7))
s13_3  <- post$rho * 13 - A_post * (1 - exp(-post$zeta * 10))

cat("\n---- Identifiable combinations ----\n")
cat("s(tau=10, Delta=+3d)  2.5% / 50% / 97.5%:\n")
print(quantile(s10_3, c(0.025, 0.5, 0.975)))
cat("s(tau=13, Delta=+3d)  2.5% / 50% / 97.5%:\n")
print(quantile(s13_3, c(0.025, 0.5, 0.975)))

cat("\nDone.\n")
