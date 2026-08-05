# fit_eip_hill.R
# Driver for eip_fit_hill.stan -- Hill-decay variant of the EIP slowdown
# fit. Mirrors fit_eip.R structure, sources priors from
# priors_atn_main.R, and saves the stanfit + JSON prior record + standard
# diagnostic plots.

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
# When run via `Rscript fit_eip_hill.R`, change the working directory to
# the script's location so all relative paths resolve. When sourced
# interactively (source("fit_eip_hill.R") in RStudio / a live R session),
# leave the working directory alone -- the user is presumably already in
# the right place, and silently moving them elsewhere is surprising.
script_dir <- tryCatch({
  args <- commandArgs(trailingOnly = FALSE)
  fn   <- sub("--file=", "", args[grep("--file=", args)])
  if (length(fn) >= 1 && nzchar(fn[1])) {
    candidate <- dirname(normalizePath(fn[1], mustWork = FALSE))
    if (dir.exists(candidate)) candidate else NULL
  } else {
    NULL
  }
}, error = function(e) NULL)

if (!is.null(script_dir) && script_dir != getwd()) {
  message("Setting working directory to script location: ", script_dir)
  tryCatch(setwd(script_dir),
           error = function(e) {
             warning("Could not setwd to '", script_dir,
                     "'; staying in '", getwd(), "'")
           })
}

data_path  <- "./exp_data/SPZ_pooled_summary_APR26.csv"
model_path <- "eip_fit_hill.stan"

# ---- Priors (single source of truth) ----
source("priors_atn_main.R")
prior_data <- prior_stan_data(priors)

# ---- Gauss-Legendre nodes/weights for the Hill J(x;n) integral ----
# Computed once via the Golub-Welsch algorithm (eigendecomposition of the
# Jacobi tridiagonal matrix), no external package required. Mapped from
# the standard GL interval [-1, 1] to [0, 1] so the weights sum to 1
# and integrate(f, 0, 1) ~ sum(gl_weights_01 * f(gl_nodes_01)).
#
# Replaces the inc_beta closed form in eip_fit_hill.stan, whose
# autodiff cost (continued-fraction expansion + digamma gradients) was
# the per-leapfrog bottleneck. M = 32 gives ~12 digits of accuracy on
# the (x, n) range we sample, far below any noise floor.
gauss_legendre_01 <- function(M) {
  i        <- seq_len(M - 1)
  off_diag <- i / sqrt(4 * i^2 - 1)
  Tmat     <- matrix(0, M, M)
  Tmat[cbind(i,     i + 1)] <- off_diag
  Tmat[cbind(i + 1, i    )] <- off_diag
  eig <- eigen(Tmat, symmetric = TRUE)
  ord <- order(eig$values)
  list(
    gl_nodes_01   = (eig$values[ord] + 1) / 2,
    gl_weights_01 = (2 * eig$vectors[1, ord]^2) / 2
  )
}
gl <- gauss_legendre_01(32L)
stopifnot(all.equal(sum(gl$gl_weights_01), 1, tolerance = 1e-12))

# ---- Data ----
dat <- read.csv(data_path)

stopifnot(
  all(c("exp_id", "post_time", "day",
        "n_CTL", "pos_CTL", "n_200mg", "pos_200mg") %in% names(dat)),
  all(dat$pos_CTL   <= dat$n_CTL),
  all(dat$pos_200mg <= dat$n_200mg)
  # No sign constraint on post_time: negative values mean exposure
  # before infection and are handled by the Delta < 0 branch in
  # eip_fit_hill.stan's likelihood.
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
    Delta_r = 10L,
    # Gauss-Legendre nodes / weights on [0, 1] for hill_J inside Stan.
    gl_M          = length(gl$gl_nodes_01),
    gl_nodes_01   = gl$gl_nodes_01,
    gl_weights_01 = gl$gl_weights_01
  ),
  prior_data
)

# ---- Compile ----
mod <- stan_model(file = model_path)

# ---- Initial values ----
# Stan's default U(-2, 2) initialisation on the unconstrained scale puts
# log_nH anywhere in (0.135, 7.39), i.e. nH up to ~1620, where
# inc_beta(1/n, 1-1/n, z) is pathologically slow (one argument near 0,
# the other near 1; continued fraction barely converges). The result
# was the chain spending all its time in the first NUTS transition --
# 24 hr without progress past iteration 1.
#
# Explicit init at the prior modes keeps the warmup well-conditioned,
# with small per-chain perturbations so the chains don't all start at
# exactly the same point (helps R-hat diagnose convergence).
n_exp <- length(unique(dat$exp_idx))
init_fn <- function() {
  list(
    log_eip_excess = priors$log_eip_excess$mean +
                     rnorm(1, 0, 0.10),
    r0_frac        = max(0.05, min(0.95,
                          priors$r0_frac_eip$alpha /
                          (priors$r0_frac_eip$alpha +
                           priors$r0_frac_eip$beta) +
                          runif(1, -0.05, 0.05))),
    log_s_half     = priors$log_s_half_eip$mean +
                     rnorm(1, 0, 0.10),
    log_nH         = priors$log_nH_eip$mean +
                     rnorm(1, 0, 0.10),
    sigma_exp      = max(0.05, 0.3 + runif(1, -0.05, 0.05)),
    z_exp          = rnorm(n_exp, 0, 0.3)
  )
}

# ---- Sample ----
fit <- sampling(
  mod,
  data    = stan_data,
  seed    = 2025,
  chains  = 4,
  iter    = 3000,
  warmup  = 1000,
  init    = init_fn,
  control = list(adapt_delta = 0.95, max_treedepth = 12),
  refresh = 500
)

# ---- Diagnostics ----
cat("\n---- HMC diagnostics ----\n")
check_hmc_diagnostics(fit)

cat("\n---- Key parameter summary ----\n")
key_vars <- c("log_eip_excess", "eip_excess", "eip_baseline",
              "rho", "rho0", "r0_frac", "B",
              "log_s_half", "s_half",
              "log_nH", "nH",
              "sigma_exp")
print(summary(fit, pars = key_vars)$summary)

cat("\n---- OLRE offsets (beta_exp) ----\n")
print(summary(fit, pars = "beta_exp")$summary)

# ---- Pairs plot ----
# Use the unconstrained primary parameters so the geometry HMC sees is visible.
post_arr <- as.array(fit, pars = c("log_eip_excess", "r0_frac",
                                   "log_s_half", "log_nH", "sigma_exp"))
p_pairs <- mcmc_pairs(post_arr,
                      off_diag_args = list(size = 0.4, alpha = 0.25))
print(p_pairs)
ggsave("eip_hill_pairs.png", p_pairs, width = 9, height = 9, dpi = 200)

# ---- Trace plot ----
p_trace <- mcmc_trace(post_arr)
ggsave("eip_hill_trace.png", p_trace, width = 10, height = 6, dpi = 200)

# ---- Posterior predictive checks ----
pp_ctl <- rstan::extract(fit, pars = "pos_ctl_rep")$pos_ctl_rep
pp_trt <- rstan::extract(fit, pars = "pos_trt_rep")$pos_trt_rep

ppc_ctl <- ppc_intervals(
  y    = dat$pos_CTL,
  yrep = pp_ctl,
  x    = seq_len(nrow(dat))
) +
  labs(title    = "Hill-decay control arm: posterior predictive",
       subtitle = "x = row index (= exp_id after re-indexing)",
       x = "row", y = "sporozoite-positive count")
ggsave("ppc_hill_ctl.png", ppc_ctl, width = 9, height = 5, dpi = 200)

ppc_trt <- ppc_intervals(
  y    = dat$pos_200mg,
  yrep = pp_trt,
  x    = seq_len(nrow(dat))
) +
  labs(title    = "Hill-decay treatment (200 mg) arm: posterior predictive",
       subtitle = "x = row index (= exp_id after re-indexing)",
       x = "row", y = "sporozoite-positive count")
ggsave("ppc_hill_trt.png", ppc_trt, width = 9, height = 5, dpi = 200)

# ---- Save fit ----
# Write to eip_fit_hill_result.rds (NOT eip_fit_hill.rds): the latter is
# the auto_write cache slot for eip_fit_hill.stan and overwriting it
# with a sampled fit produces a cryptic "no method for coercing this S4
# class to a vector" error on the next stan_model() call. Same trap as
# in fit_eip.R.
saveRDS(fit, file = "eip_fit_hill_result.rds")
write_priors_json("eip_fit_hill_priors.json", priors)
cat("\nSaved: eip_fit_hill_result.rds, eip_fit_hill_priors.json,",
    "eip_hill_pairs.png, eip_hill_trace.png, ppc_hill_ctl.png, ppc_hill_trt.png\n")

# ---- Identifiability check ----
post <- rstan::extract(fit, pars = c("rho", "rho0", "s_half", "nH", "r0_frac"))

cat("\n---- Joint posterior correlations ----\n")
cor_mat <- cor(cbind(rho     = post$rho,
                     rho0    = post$rho0,
                     s_half  = post$s_half,
                     nH      = post$nH,
                     r0_frac = post$r0_frac))
print(cor_mat)

# Effective cumulative rate at readout for Delta = +3 d, tau_obs = 10 / 13.
hill_J <- function(x, n) {
  out <- numeric(length(x))
  pos <- x > 0
  if (any(pos)) {
    z   <- x[pos]^n / (1 + x[pos]^n)
    out[pos] <- pi / (n * sin(pi / n)) * pbeta(z, 1/n, 1 - 1/n)
  }
  out
}

s_at <- function(rho, rho0, s_half, n, tau_obs, delta_val) {
  B <- rho - rho0
  out <- rho * tau_obs
  if (tau_obs > delta_val) {
    out <- out - B * s_half *
      mapply(function(rh, rh0, sh, nn) hill_J((tau_obs - delta_val)/sh, nn),
             rho, rho0, s_half, n)
  }
  out
}

s10_3 <- s_at(post$rho, post$rho0, post$s_half, post$nH, 10, 3)
s13_3 <- s_at(post$rho, post$rho0, post$s_half, post$nH, 13, 3)

cat("\n---- Identifiable combinations ----\n")
cat("s(tau=10, Delta=+3d)  2.5% / 50% / 97.5%:\n")
print(quantile(s10_3, c(0.025, 0.5, 0.975)))
cat("s(tau=13, Delta=+3d)  2.5% / 50% / 97.5%:\n")
print(quantile(s13_3, c(0.025, 0.5, 0.975)))

cat("\nDone.\n")
