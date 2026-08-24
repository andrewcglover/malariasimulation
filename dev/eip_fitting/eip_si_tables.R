# eip_si_tables.R
# Regenerates every number in the EIP-fitting tables of the ATN supplement
# (dev/reference/10AUG26_ATN_Overleaf.tex, section "EIP lengthening").
#
# Produces, printed to the console as LaTeX-ready rows:
#   Table 1  prior and posterior summaries of every fitted parameter, both forms
#   Tables 2 and 3  mean EIP duration, extension and % increase across exposure
#                   times s_0 = seq(-12, 12, 2), exponential and Hill forms
#   Tables 4 and 5  timing of exposure required to extend the EIP by at least a
#                   given amount, with the probability that it is attainable
#
# Inputs (read-only; NOT in this repository, paths are absolute by necessity):
#   Hill posterior         malariasimple_ATNs/dev/eip_fit_hill_result.rds
#   exponential posterior  malariasimple_ATNs_190526backup/dev/eip_fit_result.rds
# The Hill fit is the one consumed by malariasimulation, and is duplicated at
# dev/atn_params/eip_fit_hill_result.rds. The exponential fit exists only in the
# backup repository, its driver and Stan model having been archived here as
# fit_eip.R and eip_fit.stan.
#
# Priors are hard-coded below to match priors_atn_main.R in this folder. They are
# NOT read from that file, because prior summaries for derived quantities such as
# rho_0 and rho_min require sampling rather than the hyperparameters alone. If
# priors_atn_main.R is ever changed, change them here too.
#
# Dependencies: rstan (to read the stanfit objects), posterior (effective sample
# size for the Monte Carlo error on the attainability probabilities).
#
# Verified on 2026-08-16 to reproduce every value in tables 2 to 5 of the .tex
# exactly, and the posterior columns of table 1. Prior summaries are sampled, so
# the final digit may move between runs; see the note on N_PRIOR_TAB1 below.

suppressPackageStartupMessages({library(rstan); library(posterior)})
set.seed(7)

R_STAGES  <- 10     # number of latently-infected (EIP) compartments, fixed in the fit
EIP_FLOOR <- 7      # biological minimum on the baseline EIP, days
# Prior draws. The grid-based tables solve for the lengthened EIP at every draw,
# so their cost scales with N_PRIOR and 8000 is used, matching the posterior. The
# table 1 summaries are quantiles of the priors alone, needing no solve, so a much
# larger sample is used there: at 8000 draws the upper tail of the baseline-EIP
# prior is visibly noisy (its 97.5th centile lands anywhere between about 14.9 and
# 15.4 days against an analytic 14.91).
N_PRIOR      <- 8000
N_PRIOR_TAB1 <- 4e6

hill_rds <- "C:/Users/ag4218/Local/GitHub/malariasimple_ATNs/dev/eip_fit_hill_result.rds"
exp_rds  <- "C:/Users/ag4218/Local/GitHub/malariasimple_ATNs_190526backup/dev/eip_fit_result.rds"

# ---- suppression kernels ---------------------------------------------------
# d(s) is the cumulative progression ("distance") travelled s days after
# antimalarial exposure, i.e. the integral of the suppressed progression rate.
# Hill:        sigma(s) = s_half^eta / (s_half^eta + s^eta)
# exponential: sigma(s) = exp(-zeta * s)

# J(x; n) = int_0^x du / (1 + u^n), via the regularised incomplete beta function.
hill_J <- function(x, n) {
  out <- numeric(length(x)); k <- x > 0
  out[k] <- (pi / (n[k] * sin(pi / n[k]))) *
    pbeta(x[k]^n[k] / (1 + x[k]^n[k]), 1 / n[k], 1 - 1 / n[k])
  out
}
d_hill <- function(s, p) p$rho0 * s - (p$rho0 - p$rmin) * p$sh * hill_J(s / p$sh, p$eta)
d_exp  <- function(s, p) p$rho0 * s - ((p$rho0 - p$rmin) / p$zeta) * (1 - exp(-p$zeta * s))

# ---- lengthened EIP --------------------------------------------------------
# Mean EIP duration when exposure occurs s0 days relative to infection (s0 < 0
# meaning exposure precedes infection), found by bisection on the requirement
# that the total progression equals rho_0 * r_0. Vectorised over draws; s0 is a
# scalar. Exposure occurring after the EIP has completed (s0 >= r_0) has no
# effect, hence the final ifelse.
eip_at <- function(s0, p, d) {
  g <- if (s0 <= 0)
    function(r) d(abs(s0) + r, p) - d(abs(s0), p) - p$rho0 * p$r0
  else
    function(r) d(pmax(r - s0, 0), p) - p$rho0 * (p$r0 - s0)
  lo <- rep(1e-9, length(p$r0)); hi <- p$r0 + 60
  for (i in 1:45) {                       # 45 halvings of a ~70-day bracket
    m <- (lo + hi) / 2; neg <- g(m) < 0
    lo[neg] <- m[neg]; hi[!neg] <- m[!neg]
  }
  ifelse(s0 >= p$r0, p$r0, (lo + hi) / 2)
}
ext_at <- function(s0, p, d) eip_at(s0, p, d) - p$r0

# ---- draws -----------------------------------------------------------------
av <- function(z) as.vector(z)   # rstan returns 1-d arrays, which will not
                                 # broadcast against a matrix later on

fit_h <- readRDS(hill_rds)
ph <- rstan::extract(fit_h, pars = c("eip_baseline", "rho", "rho0", "s_half", "nH"))
post_H <- list(r0 = av(ph$eip_baseline), rho0 = av(ph$rho), rmin = av(ph$rho0),
               sh = av(ph$s_half), eta = av(ph$nH))

fit_e <- readRDS(exp_rds)
pe <- rstan::extract(fit_e, pars = c("eip_baseline", "rho", "rho0", "zeta"))
post_E <- list(r0 = av(pe$eip_baseline), rho0 = av(pe$rho), rmin = av(pe$rho0),
               zeta = av(pe$zeta))

# Priors, matching priors_atn_main.R. The Hill exponent is truncated at
# exp(0.5) in eip_fit_hill.stan, so the prior is resampled below that bound.
r0_pri   <- EIP_FLOOR + rlnorm(N_PRIOR, log(3) - 0.5, 0.75)
rho0_pri <- R_STAGES / r0_pri
eta_pri  <- rlnorm(N_PRIOR, log(5), 0.4)
while (any(eta_pri < exp(0.5))) {
  b <- eta_pri < exp(0.5); eta_pri[b] <- rlnorm(sum(b), log(5), 0.4)
}
pri_H <- list(r0 = r0_pri, rho0 = rho0_pri, rmin = rbeta(N_PRIOR, 2, 5) * rho0_pri,
              sh = rlnorm(N_PRIOR, log(1.5), 0.4), eta = eta_pri)
pri_E <- list(r0 = r0_pri, rho0 = rho0_pri, rmin = rbeta(N_PRIOR, 0.5, 0.5) * rho0_pri,
              zeta = log(2) / rlnorm(N_PRIOR, log(1.5), 0.75))

ci <- function(v, dp = 2) {
  q <- quantile(v, c(0.025, 0.5, 0.975), na.rm = TRUE)
  sprintf(paste0("%.", dp, "f [%.", dp, "f, %.", dp, "f]"), q[2], q[1], q[3])
}

# ---- Table 1: parameter summaries -----------------------------------------
# Drawn separately, and far more densely, than the priors used for tables 2 to 5;
# see the note on N_PRIOR_TAB1 above.
r0_big   <- EIP_FLOOR + rlnorm(N_PRIOR_TAB1, log(3) - 0.5, 0.75)
rho0_big <- R_STAGES / r0_big
eta_big  <- rlnorm(N_PRIOR_TAB1, log(5), 0.4)
eta_big  <- eta_big[eta_big >= exp(0.5)]
pri_H_big <- list(r0 = r0_big, rho0 = rho0_big,
                  rmin = rbeta(N_PRIOR_TAB1, 2, 5) * rho0_big,
                  sh = rlnorm(N_PRIOR_TAB1, log(1.5), 0.4), eta = eta_big)
pri_E_big <- list(r0 = r0_big, rho0 = rho0_big,
                  rmin = rbeta(N_PRIOR_TAB1, 0.5, 0.5) * rho0_big,
                  zeta = log(2) / rlnorm(N_PRIOR_TAB1, log(1.5), 0.75))

cat("\n######## TABLE 1: parameter priors and posteriors ########\n")
cat("shared      r_0 (days)        prior", ci(pri_H_big$r0),
    " exp", ci(post_E$r0), " hill", ci(post_H$r0), "\n")
cat("shared      rho_0 (/day)      prior", ci(pri_H_big$rho0),
    " exp", ci(post_E$rho0), " hill", ci(post_H$rho0), "\n")
cat("shared      sigma_beta        prior", ci(abs(rnorm(N_PRIOR_TAB1))),
    " exp", ci(av(rstan::extract(fit_e, "sigma_exp")$sigma_exp)),
    " hill", ci(av(rstan::extract(fit_h, "sigma_exp")$sigma_exp)), "\n")
cat("exponential ratio             prior", ci(pri_E_big$rmin / pri_E_big$rho0),
    " post", ci(post_E$rmin / post_E$rho0), "\n")
cat("exponential rho_min (/day)    prior", ci(pri_E_big$rmin),
    " post", ci(post_E$rmin), "\n")
cat("exponential zeta (/day)       prior", ci(pri_E_big$zeta),
    " post", ci(post_E$zeta), "\n")
cat("exponential half-life (days)  prior", ci(log(2) / pri_E_big$zeta),
    " post", ci(log(2) / post_E$zeta), "\n")
cat("hill        ratio             prior", ci(pri_H_big$rmin / pri_H_big$rho0),
    " post", ci(post_H$rmin / post_H$rho0), "\n")
cat("hill        rho_min (/day)    prior", ci(pri_H_big$rmin),
    " post", ci(post_H$rmin), "\n")
cat("hill        s_half (days)     prior", ci(pri_H_big$sh), " post", ci(post_H$sh), "\n")
cat("hill        eta               prior", ci(pri_H_big$eta), " post", ci(post_H$eta), "\n")

# ---- Tables 2 and 3: duration, extension and % increase by exposure time ---
duration_rows <- function(pri, post, d, label) {
  cat("\n######## ", label, ": duration / extension / % increase ########\n")
  for (s0 in seq(-12, 12, 2)) {
    a <- eip_at(s0, pri, d); b <- eip_at(s0, post, d)
    cat(sprintf("%d & %s & %s & %s & %s & %s & %s \\\\\n", s0,
                ci(a, 1), ci(b, 1),
                ci(a - pri$r0, 1), ci(b - post$r0, 1),
                ci(100 * (a - pri$r0) / pri$r0, 0),
                ci(100 * (b - post$r0) / post$r0, 0)))
  }
  cat("baseline r_0: prior ", ci(pri$r0), " posterior ", ci(post$r0), "\n")
}
duration_rows(pri_E, post_E, d_exp,  "TABLE 2 (exponential)")
duration_rows(pri_H, post_H, d_hill, "TABLE 3 (Hill)")

# ---- Tables 4 and 5: exposure timing required for a given extension --------
# The extension is largest at s0 = 0 and falls away on both sides, so for each
# threshold there is one crossing before infection and one after. Rather than
# invert twice per threshold, the extension is evaluated once on a grid and the
# crossings are interpolated, which is far cheaper.
GRID_PRE  <- seq(-25, 0, by = 0.10)
GRID_POST <- seq(0, 11, by = 0.05)
MIN_DRAWS <- 100   # below this the conditional summaries are too thin to report

ext_grid <- function(p, d, g) sapply(g, function(s0) ext_at(s0, p, d))

crossing <- function(M, g, target, side) {
  target <- as.vector(target); ge <- M >= target; n <- nrow(M)
  if (side == "pre") {
    k  <- max.col(ge, ties.method = "first"); lo <- pmax(k - 1, 1)
    y1 <- M[cbind(1:n, lo)]; y2 <- M[cbind(1:n, k)]
    f  <- ifelse(y2 > y1, (target - y1) / (y2 - y1), 0)
    out <- g[lo] + f * (g[k] - g[lo]); out[k == 1] <- g[1]
  } else {
    k  <- ncol(M) - max.col(ge[, ncol(M):1, drop = FALSE], ties.method = "first") + 1
    hi <- pmin(k + 1, ncol(M))
    y1 <- M[cbind(1:n, k)]; y2 <- M[cbind(1:n, hi)]
    f  <- ifelse(y1 > y2, (y1 - target) / (y1 - y2), 0)
    out <- g[k] + f * (g[hi] - g[k])
  }
  out
}

window_rows <- function(p, d, label) {
  Mpre <- ext_grid(p, d, GRID_PRE); Mpost <- ext_grid(p, d, GRID_POST)
  emax <- ext_at(0, p, d); N <- length(emax)
  cat("\n######## ", label, ": exposure timing thresholds ########\n")
  cat("maximum extension", ci(emax), " maximum %", ci(100 * emax / p$r0, 0), "\n")
  for (set in list(list(unit = "days", v = seq(0.5, 5, 0.5)),
                   list(unit = "%",    v = seq(5, 50, 5)))) {
    for (t in set$v) {
      target <- if (set$unit == "days") rep(t, N) else (t / 100) * p$r0
      ok <- emax >= target; prob <- mean(ok)
      ess <- tryCatch(posterior::ess_basic(as.numeric(ok)), error = function(e) N)
      mcse <- sqrt(prob * (1 - prob) / max(ess, 1))
      a <- crossing(Mpre, GRID_PRE, target, "pre")[ok]
      b <- crossing(Mpost, GRID_POST, target, "post")[ok]
      show <- sum(ok) >= MIN_DRAWS
      cat(sprintf("$\\geq$ %s %-4s & %.0f & %s & %s & %s \\\\  %% p=%.4f mcse=%.4f n=%d\n",
                  t, set$unit, 100 * prob,
                  if (show) ci(a) else "---",
                  if (show) ci(b) else "---",
                  if (show) ci(b - a) else "---",
                  prob, mcse, sum(ok)))
    }
  }
}
window_rows(post_E, d_exp,  "TABLE 4 (exponential)")
window_rows(post_H, d_hill, "TABLE 5 (Hill)")
