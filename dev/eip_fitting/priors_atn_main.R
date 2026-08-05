# priors_atn_main.R
# Central source of truth for the prior hyperparameters of all ATN
# main-analysis Stan models -- uni-directional, bidir shared-nH and
# bidir split-nH variants, both linear (TBA) and NB rate-thinning (TRA)
# likelihoods.  This file is sourced by:
#
#   - fit_atn_main*.R drivers, which call prior_stan_data() to build the
#     prior-hyperparameter portion of the Stan data list, and call
#     write_priors_json() to save a record alongside the fit .rds.
#   - convert_tra_to_field_tba*.R scripts, which use
#     sample_kernel_prior() to build the prior pushforward of b(s) for
#     the dashed prior 95% band.
#   - plot_atn_main*.R scripts, which do the same.
#
# Stan-side requirement: every int_fit_*_hill*.stan model file declares
# these hyperparameters as data-block inputs (NOT as transformed-data
# constants), so they receive whatever values are passed in by the
# driver.  This is the only mechanism that keeps the Stan models and the
# plotting code in sync.  Do NOT redefine these values anywhere else.
#
# Prior philosophy (2026-04-29): minimally informative on the kernel
# parameters, letting the data drive the posterior.
#   r_min       ~ Beta(0.5, 0.5)        (Jeffreys reference prior)
#   log_s_half  ~ N(log(2), 1.0)        (median 2 d, ~95% in [0.28, 14.2])
#   log_nH      ~ N(log(5), 1.0)        (median 5,   ~95% in [0.71, 35.5])
# Likelihood-side hyperparameters retained at their previous defaults:
#   mu_L        ~ N(1.0, 1.5)           (NB rate-thinning, log-rate)
#   mu_p        ~ N(2.0, 1.5)           (linear / binomial, logit)
#   sigma_exp   ~ N+(0, 1)              (OLRE scale)
#   reciprocal_phi ~ N+(0, 1)           (NB dispersion, parameterised
#                                        on 1/phi)
#
# EIP-fit hyperparameters (used by eip_fit.stan / fit_eip.R).  Naming
# convention: "excess" denotes the part of the *baseline* (i.e. control)
# mean EIP that sits above the biological floor; "extension" is reserved
# for the ATN-induced lengthening of the EIP in the exposed cohort, which
# is governed by r0_frac and log_zeta below.
#   eip_floor       = 7 days                  (hard biological minimum on
#                                              baseline EIP; lower bound
#                                              below which Pf-in-An. gambiae
#                                              transmission is essentially
#                                              never observed)
#   log_eip_excess  ~ N(log(3) - 0.5, 1.0)    (lognormal on the baseline-EIP
#                                              excess above the floor; prior
#                                              mean excess = 3 d -> prior
#                                              mean baseline EIP = 10 d;
#                                              ~95% prior interval on excess
#                                              is [0.26, 13] d, baseline
#                                              EIP [7.3, 20] d)
#   r0_frac         ~ Beta(2, 2)              (rho_0 / rho in (0, 1); the
#                                              ATN-induced fractional
#                                              slowdown of the EIP rate
#                                              immediately after exposure)
#   log_zeta        ~ N(log(log 2 / 1.5), 0.2) (recovery rate of the EIP
#                                              progression rate back to
#                                              baseline after ATN exposure;
#                                              centred on a 1.5-day half
#                                              life)
# Edit r0_frac and log_zeta here to change the *prior* on the ATN-induced
# EIP extension; the change propagates to (i) the Stan likelihood via
# prior_stan_data() and the data block of eip_fit.stan, and (ii) the
# dashed prior pushforward overlay in eip_posterior_plots.R.

priors <- list(
  r_min           = list(alpha = 0.5, beta = 0.5),
  log_s_half      = list(mean  = log(2.0), sd = 1.0),
  log_nH          = list(mean  = log(5.0), sd = 1.0),
  mu_L            = list(mean  = 1.0, sd = 1.5),
  mu_p            = list(mean  = 2.0, sd = 1.5),
  sigma_exp       = list(sd    = 1.0),
  reciprocal_phi  = list(sd    = 1.0),
  eip_floor       = 7.0,
  log_eip_excess  = list(mean  = log(3.0) - 0.5, sd = 0.75),
  r0_frac         = list(alpha = 0.5, beta = 0.5),
  log_zeta        = list(mean  = log(log(2.0) / 1.5), sd = 0.75),
  # Hill-EIP-specific regularising priors. Used ONLY by eip_fit_hill.stan
  # (the exponential model and all blocking models keep using r0_frac,
  # log_s_half, log_nH above). Tighter than their shared counterparts to
  # stabilise HMC, which suffered from extreme-tail excursions into
  # numerical-pathology regions on Hill EIP fits:
  #   log_s_half_eip ~ N(log 1.5, 0.4)  -> 95% on [0.69, 3.3] d
  #     centres on the original 1.5-d half-decay timescale of the
  #     exponential model's log_zeta, with a tighter spread that keeps
  #     the chain off implausibly long timescales.
  #   log_nH_eip     ~ N(log 5,   0.4)  -> 95% on [2.27, 11.0]
  #     stays clearly sigmoidal, well away from the n_H = 1 singularity
  #     of the J(x;n) closed form (and its associated autodiff blow-ups)
  #     and away from the step-function regime n_H >> 10.
  #   r0_frac_eip    ~ Beta(2, 5)       -> mean 0.286, 95% on [0.06, 0.61]
  #     pushes prior mass away from r0_frac -> 1, where B = rho - rho_0
  #     vanishes and the entire Hill kernel becomes invisible to the
  #     likelihood (the dominant pathology under Beta(0.5, 0.5)).
  log_s_half_eip  = list(mean  = log(1.5), sd = 0.4),
  log_nH_eip      = list(mean  = log(5.0), sd = 0.4),
  r0_frac_eip     = list(alpha = 2.0, beta = 5.0)
)
# Notes on the EIP-fit prior choices:
#   - log_eip_excess sd tightened to 0.75 (95% prior interval on
#     baseline EIP roughly [7.4, 17] d, still admitting all biologically
#     plausible Pf-in-An. gambiae EIPs).
#   - r0_frac switched to a Jeffreys Beta(0.5, 0.5) reference prior so
#     the data alone drive the post-exposure fractional rate; this is
#     the same choice as for r_min in the blocking fit.
#   - log_zeta sd widened to 0.75 (recovery half-life 95% prior interval
#     widens from ~[1.0, 2.2] d at sd 0.2 to roughly [0.34, 6.6] d at
#     sd 0.75); the half-life is now genuinely data-driven rather than
#     pinned by the prior.
#   - log_s_half and log_nH (already in the list above for the blocking
#     fit) are *reused* by eip_fit_hill.stan for the Hill-decay EIP
#     model. The blocking and Hill-EIP recovery timescales are similar
#     in magnitude (~1-2 d), so the same lognormal at log 2 with sd 1.0
#     remains minimally informative for both.

# ---- Stan data wiring -----------------------------------------------------
# Build a named list of prior hyperparameters to merge into the Stan data
# list.  Each Stan file declares the subset it actually uses; extra keys
# in the data list are silently ignored by rstan, so passing the full set
# from every driver is safe and keeps drivers symmetric.
prior_stan_data <- function(priors_ = priors) {
  list(
    r_min_alpha       = priors_$r_min$alpha,
    r_min_beta        = priors_$r_min$beta,
    log_s_half_mean   = priors_$log_s_half$mean,
    log_s_half_sd     = priors_$log_s_half$sd,
    log_nH_mean       = priors_$log_nH$mean,
    log_nH_sd         = priors_$log_nH$sd,
    mu_L_mean         = priors_$mu_L$mean,
    mu_L_sd           = priors_$mu_L$sd,
    mu_p_mean         = priors_$mu_p$mean,
    mu_p_sd           = priors_$mu_p$sd,
    sigma_exp_sd      = priors_$sigma_exp$sd,
    reciprocal_phi_sd = priors_$reciprocal_phi$sd,
    eip_floor             = priors_$eip_floor,
    log_eip_excess_mean   = priors_$log_eip_excess$mean,
    log_eip_excess_sd     = priors_$log_eip_excess$sd,
    r0_frac_alpha         = priors_$r0_frac$alpha,
    r0_frac_beta          = priors_$r0_frac$beta,
    log_zeta_mean         = priors_$log_zeta$mean,
    log_zeta_sd           = priors_$log_zeta$sd,
    # Hill-EIP-specific regularising priors (read only by eip_fit_hill.stan).
    log_s_half_eip_mean   = priors_$log_s_half_eip$mean,
    log_s_half_eip_sd     = priors_$log_s_half_eip$sd,
    log_nH_eip_mean       = priors_$log_nH_eip$mean,
    log_nH_eip_sd         = priors_$log_nH_eip$sd,
    r0_frac_eip_alpha     = priors_$r0_frac_eip$alpha,
    r0_frac_eip_beta      = priors_$r0_frac_eip$beta
  )
}

# ---- Prior pushforward helpers --------------------------------------------
# Sample (r_min, s_half, nH) directly from the priors and return a list
# of vectors usable by convert/plot scripts.  variant controls the joint
# structure of the kernel parameters:
#   "uni"          -> single (s_half, nH) per draw
#   "bidir_shared" -> (s_half_pre, s_half_post) i.i.d., nH shared
#   "bidir_split"  -> (s_half_pre, s_half_post, nH_pre, nH_post) all i.i.d.
# The shared shape b_max = 1 - r_min is included for convenience.
sample_kernel_prior <- function(n_prior, variant, priors_ = priors,
                                seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  variant <- match.arg(variant, c("uni", "bidir_shared", "bidir_split"))

  out <- list()
  out$r_min <- rbeta(n_prior, priors_$r_min$alpha, priors_$r_min$beta)
  out$b_max <- 1 - out$r_min

  rln_s <- function() exp(rnorm(n_prior, priors_$log_s_half$mean,
                                priors_$log_s_half$sd))
  rln_n <- function() exp(rnorm(n_prior, priors_$log_nH$mean,
                                priors_$log_nH$sd))

  if (variant == "uni") {
    out$s_half <- rln_s()
    out$nH     <- rln_n()
  } else if (variant == "bidir_shared") {
    out$s_half_pre  <- rln_s()
    out$s_half_post <- rln_s()
    out$nH          <- rln_n()
  } else {  # bidir_split
    out$s_half_pre  <- rln_s()
    out$s_half_post <- rln_s()
    out$nH_pre      <- rln_n()
    out$nH_post     <- rln_n()
  }
  out
}

# ---- JSON record ----------------------------------------------------------
# Write the canonical priors list to a JSON file alongside the fit object.
# Round-trips with jsonlite::fromJSON for downstream inspection / re-plotting.
write_priors_json <- function(path, priors_ = priors) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("jsonlite is required to write the priors JSON record; install it.")
  }
  jsonlite::write_json(priors_, path = path,
                       pretty = TRUE, auto_unbox = TRUE,
                       digits = 12)
  invisible(path)
}
