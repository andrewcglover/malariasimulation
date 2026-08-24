# blocking_si_tables.R
# ---------------------------------------------------------------------------
# Regenerates every number quoted in the transmission-blocking section of
# dev/reference/10AUG26_ATN_Overleaf.tex, and prints LaTeX-ready table rows.
#
# THREE blocking quantities are described in that section:
#   TRA        proportional reduction in the expected OOCYST COUNT.  Fitted to
#              mosquito-level counts (NB rate-thinning likelihood).  Assumed
#              common to laboratory and field.
#   TBA_lab    proportional reduction in the PROBABILITY of infection under the
#              saturating conditions of the assay.  Fitted to per-row prevalence
#              (binomial likelihood, p_trt = (1 - b) p_ctl).
#   TBA_field  the same probability reduction under field conditions.  NOT
#              fitted; obtained by pushing TRA through the Bompard map F().
# TRA and TBA_lab share the Hill blocking function B() and ALL of its priors, so
# their prior pushforwards are identical by construction; the TBA_field prior is
# that same pushforward mapped through F().
#
# Tables produced:
#   1. tab:blk_priors_posteriors  parameters: prior | TRA posterior | TBA_lab posterior
#   2. tab:blk_probability        s_0 | lab prior | TRA | TBA_lab | field prior | TBA_field
#   3. tab:blk_window             exposure window by threshold, three scales, POSTERIOR ONLY
#      (posterior-only because under the shared Beta(0.5,0.5) prior on b_max a
#       threshold is unattainable in most prior draws: Pr = 0.50/0.33/0.21/0.14/0.06
#       at 50/75/90/95/99%, so prior windows would condition on a small minority.)
#
# Companion to eip_si_tables.R, which does the same for the EIP tables.
#
# INPUTS (both in this repo):
#   dev/atn_params/atn_fit_tra_main_bidir_splitnH.rds   (TRA, NB rate-thinning)
#   dev/atn_params/atn_fit_tba_main_bidir_splitnH.rds   (TBA_lab, binomial)
# Priors mirror ../malariasimple_ATNs/dev/priors_atn_main.R; change both together.
#
# NOTES
#   - Prior summaries in table 1 are EXACT quantiles (all marginals are Beta,
#     lognormal, half-normal or normal), so they carry no Monte Carlo error.
#     The B(s) pushforward in tables 2 and 3 is a product of random variables
#     and IS sampled, at N_PRIOR draws with a fixed seed.
#   - rstan::extract returns 1-d arrays that will not broadcast against a
#     matrix; as.vector() everything before arithmetic.
#   - Bompard parameters are those used by convert_tra_to_field_tba_bidir_
#     splitnH.R and by malariasimulation (use_bompard = TRUE), which applies the
#     map to BOTH the pre-infection FOI reduction and the post-infection
#     clearance probability (R/biting_process.R:398 and :417).
# ---------------------------------------------------------------------------

suppressPackageStartupMessages(library(rstan))

ROOT     <- "C:/Users/ag4218/Local/GitHub/malariasimulation/dev/atn_params"
TRA_PATH <- file.path(ROOT, "atn_fit_tra_main_bidir_splitnH.rds")
LAB_PATH <- file.path(ROOT, "atn_fit_tba_main_bidir_splitnH.rds")
N_PRIOR  <- 8000L      # matches the posterior draw count
SEED     <- 2025L

# Bompard / Challenger field parameters (wild An. gambiae, Burkina Faso).
M_BOMP <- 1.57e-4
NU_BOMP <- 4.95e-6

# Priors (mirror of priors_atn_main.R) --------------------------------------
PR <- list(
  r_min      = c(alpha = 0.5, beta = 0.5),
  log_s_half = c(mean = log(2), sd = 1),
  log_nH     = c(mean = log(5), sd = 1),
  mu_L       = c(mean = 1, sd = 1.5),    # NB rate-thinning (TRA), log-rate
  mu_p       = c(mean = 2, sd = 1.5),    # binomial (TBA_lab), logit
  sigma_exp  = c(sd = 1),
  recip_phi  = c(sd = 1)
)

# ---- helpers --------------------------------------------------------------
qs   <- function(x) unname(quantile(x, c(0.025, 0.5, 0.975), na.rm = TRUE))
fmt3 <- function(v, d = 3) sprintf(paste0("%.", d, "g [%.", d, "g, %.", d, "g]"),
                                   v[2], v[1], v[3])
fmt1 <- function(v, d = 1) sprintf(paste0("%.", d, "f [%.", d, "f, %.", d, "f]"),
                                   v[2], v[1], v[3])

# The general blocking function B(s_0): two-sided Hill, shared peak b_max.
B_fun <- function(s0, p) {
  a <- abs(s0)
  if (s0 > 0) p$b_max * p$s_half_post^p$nH_post / (p$s_half_post^p$nH_post + a^p$nH_post)
  else        p$b_max * p$s_half_pre^p$nH_pre   / (p$s_half_pre^p$nH_pre   + a^p$nH_pre)
}

# The laboratory-to-field map F(): Bompard, expm1-stable.
F_map <- function(tra, m = M_BOMP, nu = NU_BOMP) {
  t <- pmin(pmax(tra, 0), 1 - 1e-12)
  log_a <- nu * (log(nu) - log(nu + m))
  log_b <- nu * (log(nu) - log(nu + m * (1 - t)))
  out   <- exp(log_a) * expm1(log_b - log_a) / -expm1(log_a)
  out[t <= 0] <- 0
  pmin(pmax(out, 0), 1)
}
F_inv <- function(z) vapply(z, function(zz) {
  if (zz <= 0) return(0); if (zz >= 1) return(1)
  uniroot(function(x) F_map(x) - zz, c(0, 1 - 1e-12), tol = 1e-12)$root
}, numeric(1))

kern_pars <- c("b_max", "s_half_pre", "s_half_post", "nH_pre", "nH_post")

# ---- posterior draws ------------------------------------------------------
fit_tra <- readRDS(TRA_PATH)
fit_lab <- readRDS(LAB_PATH)
tra <- lapply(rstan::extract(fit_tra, pars = c(kern_pars, "r_min", "mu_L",
                                               "sigma_exp", "phi")), as.vector)
lab <- lapply(rstan::extract(fit_lab, pars = c(kern_pars, "r_min", "mu_p",
                                               "sigma_exp")), as.vector)

# ---- prior draws (pushforward only; shared by TRA and TBA_lab) ------------
set.seed(SEED)
prior <- list(
  b_max       = 1 - rbeta(N_PRIOR, PR$r_min["alpha"], PR$r_min["beta"]),
  s_half_pre  = exp(rnorm(N_PRIOR, PR$log_s_half["mean"], PR$log_s_half["sd"])),
  s_half_post = exp(rnorm(N_PRIOR, PR$log_s_half["mean"], PR$log_s_half["sd"])),
  nH_pre      = exp(rnorm(N_PRIOR, PR$log_nH["mean"], PR$log_nH["sd"])),
  nH_post     = exp(rnorm(N_PRIOR, PR$log_nH["mean"], PR$log_nH["sd"]))
)

# ===========================================================================
# DIAGNOSTICS
# ===========================================================================
cat("\n===== SAMPLER DIAGNOSTICS =====\n")
diag_one <- function(fit, nm) {
  sa <- summary(fit)$summary
  # b_max / r_min can return NaN n_eff when the posterior is pinned to a boundary
  par_rows <- setdiff(rownames(sa)[!is.na(sa[, "n_eff"])], "lp__")
  sp <- get_sampler_params(fit, inc_warmup = FALSE)
  cat(sprintf("%s: draws %d | max Rhat %.4f | min ESS(params) %.0f (%s) | ESS(lp__) %.0f\n",
              nm, nrow(as.matrix(fit)), max(sa[, "Rhat"], na.rm = TRUE),
              min(sa[par_rows, "n_eff"]),
              par_rows[which.min(sa[par_rows, "n_eff"])], sa["lp__", "n_eff"]))
  cat(sprintf("   divergences %d | deepest treedepth %d | E-BFMI %s\n",
              sum(sapply(sp, function(x) sum(x[, "divergent__"]))),
              max(sapply(sp, function(x) max(x[, "treedepth__"]))),
              paste(sprintf("%.2f", rstan::get_bfmi(fit)), collapse = ", ")))
}
diag_one(fit_tra, "TRA    ")
diag_one(fit_lab, "TBA_lab")

# ===========================================================================
# TABLE 1: parameter priors and posteriors (prior | TRA | TBA_lab)
# ===========================================================================
cat("\n===== TABLE 1: blocking parameters =====\n")
p <- c(0.025, 0.5, 0.975)
pri <- list(
  b_max    = rev(1 - qbeta(p, PR$r_min["alpha"], PR$r_min["beta"])),
  s_half   = qlnorm(p, PR$log_s_half["mean"], PR$log_s_half["sd"]),
  nH       = qlnorm(p, PR$log_nH["mean"], PR$log_nH["sd"]),
  mu_L     = qnorm(p, PR$mu_L["mean"], PR$mu_L["sd"]),
  exp_mu_L = exp(qnorm(p, PR$mu_L["mean"], PR$mu_L["sd"])),
  mu_p     = qnorm(p, PR$mu_p["mean"], PR$mu_p["sd"]),
  ilg_mu_p = plogis(qnorm(p, PR$mu_p["mean"], PR$mu_p["sd"])),
  sigma    = qnorm(0.5 + p / 2, 0, PR$sigma_exp["sd"]),
  phi      = rev(1 / qnorm(0.5 + p / 2, 0, PR$recip_phi["sd"]))
)
row <- function(lbl, pr, a, b) cat(sprintf(" & %s & %s & %s & %s \\\\\n", lbl, pr, a, b))
f4 <- function(v) sprintf("%.4f [%.4f, %.4f]", v[2], v[1], v[3])
row("$b_{\\max}$, peak blocking probability", sprintf("%.3f [%.4f, %.4f]", pri$b_max[2], pri$b_max[1], pri$b_max[3]),
    f4(qs(tra$b_max)), f4(qs(lab$b_max)))
row("$s^-_{1/2}$, pre-infection half-life (days)",  fmt3(pri$s_half), fmt3(qs(tra$s_half_pre)),  fmt3(qs(lab$s_half_pre)))
row("$\\eta^-$, pre-infection Hill coefficient",    fmt3(pri$nH),     fmt3(qs(tra$nH_pre)),      fmt3(qs(lab$nH_pre)))
row("$s^+_{1/2}$, post-infection half-life (days)", fmt3(pri$s_half), fmt3(qs(tra$s_half_post)), fmt3(qs(lab$s_half_post)))
row("$\\eta^+$, post-infection Hill coefficient",   fmt3(pri$nH),     fmt3(qs(tra$nH_post)),     fmt3(qs(lab$nH_post)))
cat("\\midrule\n")
row("$\\mu_L$, control-arm log burden",         fmt3(pri$mu_L),     fmt3(qs(tra$mu_L)), "---")
row("$\\exp(\\mu_L)$, control-arm burden (oocysts)", fmt3(pri$exp_mu_L), fmt3(qs(exp(tra$mu_L))), "---")
row("$\\mu_p$, control-arm log-odds of infection", fmt3(pri$mu_p), "---", fmt3(qs(lab$mu_p)))
row("$\\mathrm{logit}^{-1}(\\mu_p)$, control-arm prevalence", fmt3(pri$ilg_mu_p), "---", fmt3(qs(plogis(lab$mu_p))))
row("$\\sigma_\\beta$", fmt3(pri$sigma), fmt3(qs(tra$sigma_exp)), fmt3(qs(lab$sigma_exp)))
row("$\\phi$, dispersion", fmt3(pri$phi), fmt3(qs(tra$phi)), "---")

# ===========================================================================
# TABLE 2: blocking probability across exposure times
# ===========================================================================
cat("\n===== TABLE 2: blocking probability by exposure time =====\n")
S0 <- seq(-6, 6, 1)
for (s0 in S0) {
  pr_lab <- B_fun(s0, prior)
  cat(sprintf("%g & %s & %s & %s & %s & %s \\\\\n", s0,
              fmt1(100 * qs(pr_lab)),
              fmt1(100 * qs(B_fun(s0, tra))),
              fmt1(100 * qs(B_fun(s0, lab))),
              fmt1(100 * qs(F_map(pr_lab))),
              fmt1(100 * qs(F_map(B_fun(s0, tra))))))
}

# ===========================================================================
# TABLE 3: exposure window by blocking threshold (posterior only, 3 scales)
# ===========================================================================
cat("\n===== TABLE 3: exposure window by threshold =====\n")
# B(s) >= theta  <=>  s <= s_half * ((b_max - theta)/theta)^(1/nH)
edge <- function(theta, b_max, s_h, nH) {
  ok <- b_max >= theta
  out <- rep(NA_real_, length(b_max))
  out[ok] <- s_h[ok] * ((b_max[ok] - theta) / theta)^(1 / nH[ok])
  out
}
thresholds <- c(0.25, 0.50, 0.75, 0.90, 0.95, 0.99)
window_block <- function(pars, scale) {
  for (th in thresholds) {
    th_b <- if (scale == "field") F_inv(th) else th
    pre  <- edge(th_b, pars$b_max, pars$s_half_pre,  pars$nH_pre)
    post <- edge(th_b, pars$b_max, pars$s_half_post, pars$nH_post)
    keep <- !is.na(pre)
    cat(sprintf("$\\geq$ %g\\%% & %s & %s & %s \\\\ %% attainable %.1f%%, B-threshold %.4f\n",
                100 * th,
                fmt1(-rev(qs(pre[keep])), 2), fmt1(qs(post[keep]), 2),
                fmt1(qs((pre + post)[keep]), 2),
                100 * mean(keep), th_b))
  }
}
cat("-- TRA --\n");       window_block(tra, "lab")
cat("-- TBA_lab --\n");   window_block(lab, "lab")
cat("-- TBA_field --\n"); window_block(tra, "field")

cat("\nPrior attainability (why table 3 is posterior-only):\n")
for (th in thresholds)
  cat(sprintf("  threshold %.2f : prior Pr(b_max >= th) = %.3f\n", th, 1 - pbeta(th, 0.5, 0.5)))

# ===========================================================================
# In-text quantities
# ===========================================================================
cat("\n===== IN-TEXT NUMBERS =====\n")
cat(" s0 |        TRA         |      TBA_lab       |     TBA_field\n")
for (s0 in c(-5, -4, -3, -2, -1, -0.5, 0, 0.5, 1, 2, 3, 5)) {
  cat(sprintf("%5g | %-18s | %-18s | %s\n", s0,
              fmt1(100 * qs(B_fun(s0, tra))),
              fmt1(100 * qs(B_fun(s0, lab))),
              fmt1(100 * qs(F_map(B_fun(s0, tra))))))
}
cat("\nF() at fixed TRA values (deterministic):\n")
for (t in c(0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 0.999, 1))
  cat(sprintf("  TRA %.3f -> TBA_field %.4f\n", t, F_map(t)))
cat(sprintf("\nField: m/nu = %.1f ; P(infected per bite) = %.3g\n",
            M_BOMP / NU_BOMP, 1 - (NU_BOMP / (NU_BOMP + M_BOMP))^NU_BOMP))
cat(sprintf("Field burden among infected mosquitoes = %.2f oocysts\n",
            M_BOMP / (1 - (NU_BOMP / (NU_BOMP + M_BOMP))^NU_BOMP)))
cat(sprintf("Halving the rate: P(infected) falls %.1f%%, burden among infected %.2f -> %.2f\n",
            100 * (1 - (1 - (NU_BOMP / (NU_BOMP + M_BOMP / 2))^NU_BOMP) /
                     (1 - (NU_BOMP / (NU_BOMP + M_BOMP))^NU_BOMP)),
            M_BOMP / (1 - (NU_BOMP / (NU_BOMP + M_BOMP))^NU_BOMP),
            (M_BOMP / 2) / (1 - (NU_BOMP / (NU_BOMP + M_BOMP / 2))^NU_BOMP)))

cat("\nHalf-life asymmetry (posterior of the difference):\n")
for (nm in c("tra", "lab")) {
  pp <- get(nm)
  d  <- pp$s_half_pre - pp$s_half_post
  cat(sprintf("  %-8s s^-_{1/2} - s^+_{1/2} = %s d ; Pr(pre > post) = %.4f\n",
              nm, fmt3(qs(d)), mean(d > 0)))
}

cat("\nAssay burdens (from the mosquito-level CSV):\n")
mos <- read.csv("C:/Users/ag4218/Local/GitHub/malariasimple_ATNs/dev/exp_data/mos_level_int_summary_pre_post_full.csv")
ctl <- mos$int[mos$group == "control"]
cat(sprintf("  control: n = %d, positive = %d (%.1f%%), mean = %.1f, mean among infected = %.1f\n",
            length(ctl), sum(ctl > 0), 100 * mean(ctl > 0), mean(ctl), mean(ctl[ctl > 0])))
