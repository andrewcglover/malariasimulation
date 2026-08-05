// eip_fit_hill.stan
// ATN EIP-slowdown model fit, Hill-decay variant of eip_fit.stan.
//
// Mechanism (for fresh nets, tau = 0 in the PDF):
//   After exposure at time Delta relative to infection, the progression
//   rate through the EIP is depressed by a Hill-shaped suppression that
//   recovers to the baseline rho with characteristic time s_half:
//
//     rho(s) = rho - (rho - rho_0) * sigma(s),
//     sigma(s) = s_half^n / (s_half^n + s^n),    s = time since exposure
//
//   sigma(0) = 1 (full impairment at exposure), sigma(infty) = 0
//   (full recovery), and for n > 1 the suppression is sigmoidal in s.
//
// Time-change argument:
//   EIP completion has the distribution of the hitting time of
//   Gamma(shape = Delta_r, rate = 1) by the integrated rate
//   s(tau_obs; Delta) (same Erlang accumulator as eip_fit.stan):
//
//     P(EIP completed by tau_obs | Delta)
//       = F_Gamma(s(tau_obs; Delta); Delta_r, 1).
//
//   Closed forms (B = rho - rho_0):
//
//   Case (i)  Delta >= 0, tau_obs > Delta:
//     s(tau_obs; Delta) = rho * tau_obs - B * s_half * J(x; n),
//                          x = (tau_obs - Delta) / s_half.
//
//   Case (ii) Delta >= 0, tau_obs <= Delta:
//     s(tau_obs; Delta) = rho * tau_obs.
//
//   Case (iii) Delta < 0 (exposure before infection):
//     s(tau_obs; Delta) = rho * tau_obs
//                       - B * s_half * (J(x_hi; n) - J(x_lo; n)),
//     with x_hi = (tau_obs - Delta) / s_half and x_lo = -Delta / s_half.
//     The J(x_lo; n) term subtracts the integrated suppression that
//     occurred during the |Delta| days between exposure and infection,
//     leaving only the post-infection portion of the integral.
//
//   J(x; n) = int_0^x du / (1 + u^n).
//
//   Closed form (mathematically): via the regularised incomplete beta
//   function, valid for n > 1:
//     J(x; n) = (pi / (n * sin(pi/n)))
//             * inc_beta(1/n, 1 - 1/n, x^n / (1 + x^n)).
//
//   Numerical implementation (used here): TWO-PIECE M-point Gauss-
//   Legendre quadrature on [0, 1] after the substitution u = x*v.
//   The interval is split at v_star = 1/x (the location of the sharp
//   transition of the integrand 1/(1+(xv)^n)), and M GL nodes are
//   applied to each subinterval. With M = 32 this gives ~12 digits of
//   accuracy uniformly in (x, n) over the range we sample. A naive
//   single-piece M = 32 GL on [0, 1] was tested and left ~1e-3 error
//   at (x=10, n=10) because few nodes landed near the transition;
//   splitting fixes this entirely.
//
//   GL nodes/weights are computed once on the R side via Golub-Welsch
//   (eigendecomposition of the Jacobi tridiagonal) and passed in via
//   the data block. Per-call cost is fixed at 2M pow() evaluations,
//   typically 5-50x faster in autodiff than the inc_beta closed
//   form, whose continued-fraction gradient with respect to a, b, x
//   dominated the original implementation's leapfrog cost (~9 sec/
//   step for N = 25 rows).
//
// Identifiability + numerical-stability constraint:
//   log_nH is constrained to >= 0.5, i.e. nH >= exp(0.5) ~ 1.65. The
//   biological motivation is that nH = 1 is the sigmoidality boundary
//   (Hill recovery becomes hyperbolic), and the numerical motivation
//   is that Stan's `inc_beta(1/n, 1 - 1/n, z)` continued-fraction
//   expansion is slow and unstable when either argument is close to 0,
//   which happens at both nH -> 1+ (b -> 0+) and nH -> infty (a -> 0+).
//   Lifting the lower bound from 0 to 0.5 keeps autodiff gradients
//   finite during warmup, even when Stan's default U(-2, 2) init
//   draws perturb the chain away from the prior mode.
//
//   With log_nH_eip ~ N(log 5, 0.4), only ~0.3% of prior mass falls
//   below log_nH = 0.5 (truncation point at z = (0.5 - log 5)/0.4
//   ~ -2.77 std devs), so the truncation is biologically negligible.
//
// Likelihood: paired Binomial per row with a shared logit-scale OLRE
// offset beta_exp[exp_id] applied identically to control and 200mg arms,
// matching eip_fit.stan.

functions {
  // J(x; n) = int_0^x du / (1 + u^n) via TWO-PIECE M-point Gauss-Legendre
  // quadrature on [0, 1] (after the substitution u = x*v). The GL nodes
  // and weights are passed in pre-mapped to [0, 1].
  //
  // The integrand 1/(1 + (xv)^n) has a sharp transition at v = 1/x for
  // large x, n. Single-piece M = 32 GL on [0, 1] left ~1e-3 error at
  // (x=10, n=10) because few nodes lie near the transition. Splitting
  // [0, 1] at v_star = 1/x (or v_star = 0.5 if x <= 1) puts the
  // transition exactly at a subinterval boundary, giving each side its
  // own dense GL grid and dropping the worst-case error to ~1e-12
  // across the entire (x, n) range we sample. Cost is fixed at 2M
  // pow() evaluations per call.
  real hill_J(real x, real n,
              vector gl_nodes_01, vector gl_weights_01) {
    if (x <= 0) return 0.0;
    int M = num_elements(gl_nodes_01);

    // Split point: at the transition for x > 1, midpoint otherwise.
    real v_star = (x > 1.0) ? (1.0 / x) : 0.5;
    real width_hi = 1.0 - v_star;

    // Lower piece: integrate over [0, v_star].
    real total_lo = 0.0;
    for (i in 1:M) {
      real v = v_star * gl_nodes_01[i];
      total_lo += gl_weights_01[i] / (1.0 + pow(x * v, n));
    }
    total_lo *= v_star;

    // Upper piece: integrate over [v_star, 1].
    real total_hi = 0.0;
    for (i in 1:M) {
      real v = v_star + width_hi * gl_nodes_01[i];
      total_hi += gl_weights_01[i] / (1.0 + pow(x * v, n));
    }
    total_hi *= width_hi;

    return x * (total_lo + total_hi);
  }
}

data {
  int<lower=1> N;                              // rows
  int<lower=1> K;                              // unique experiments
  array[N] int<lower=1, upper=K> exp_id;       // 1..K
  vector[N] delta;                             // days: post_time / 24,
                                               // negative = exposure before infection
  vector<lower=0>[N] tau_obs;                  // days post infection (10 or 13)
  array[N] int<lower=0> n_ctl;
  array[N] int<lower=0> pos_ctl;
  array[N] int<lower=0> n_trt;
  array[N] int<lower=0> pos_trt;
  int<lower=1> Delta_r;                        // Erlang shape (= 10)

  // EIP prior hyperparameters (sourced from priors_atn_main.R).  The
  // baseline-EIP-floor reparameterisation matches eip_fit.stan; the
  // Hill kernel parameters use the *_eip suffixed regularising priors
  // dedicated to this model (NOT the shared log_s_half / log_nH used
  // by the blocking fits, and NOT the r0_frac used by the exponential
  // EIP fit). This keeps the exp pipeline's priors entirely untouched.
  real<lower=0> eip_floor;
  real log_eip_excess_mean;
  real<lower=0> log_eip_excess_sd;
  real<lower=0> r0_frac_eip_alpha;
  real<lower=0> r0_frac_eip_beta;
  real log_s_half_eip_mean;
  real<lower=0> log_s_half_eip_sd;
  real log_nH_eip_mean;
  real<lower=0> log_nH_eip_sd;

  // Gauss-Legendre nodes and weights on [0, 1] for the hill_J integral.
  // Pre-computed in fit_eip_hill.R via Golub-Welsch and passed in.
  // M = 32 gives ~12 digits of accuracy on the integrand range we sample.
  int<lower=1> gl_M;
  vector<lower=0, upper=1>[gl_M] gl_nodes_01;
  vector<lower=0>[gl_M] gl_weights_01;
}

transformed data {
  real shape_r = Delta_r;
}

parameters {
  real log_eip_excess;                 // log of baseline-EIP excess above floor (days)
  real<lower=0, upper=1> r0_frac;      // rho_0 / rho in (0, 1); ATN-induced
                                       // fractional slowdown at exposure
  real log_s_half;                     // log of Hill half-recovery time (days)
  real<lower=0.5> log_nH;              // log of Hill exponent; nH = exp(log_nH) >= ~1.65
                                       // (lifted from 0 to 0.5 for numerical stability;
                                       //  see header for rationale)
  real<lower=0> sigma_exp;             // SD of OLRE on logit scale
  vector[K] z_exp;                     // non-centred OLRE draws
}

transformed parameters {
  real<lower=0> eip_excess   = exp(log_eip_excess);
  real<lower=0> eip_baseline = eip_floor + eip_excess;
  real<lower=0> rho          = shape_r / eip_baseline;
  real<lower=0> s_half       = exp(log_s_half);
  real<lower=1.6> nH         = exp(log_nH);   // = exp(>=0.5) >= ~1.65
  real<lower=0> rho0         = r0_frac * rho;
  real<lower=0> B            = rho - rho0;
  vector[K] beta_exp         = sigma_exp * z_exp;
}

model {
  // Priors (all hyperparameters from priors_atn_main.R via the data block).
  // Note the *_eip suffixes on the kernel-parameter priors -- these are
  // the regularised Hill-EIP-only priors, NOT the shared blocking priors.
  log_eip_excess ~ normal(log_eip_excess_mean, log_eip_excess_sd);
  r0_frac        ~ beta(r0_frac_eip_alpha, r0_frac_eip_beta);
  log_s_half     ~ normal(log_s_half_eip_mean, log_s_half_eip_sd);
  log_nH         ~ normal(log_nH_eip_mean, log_nH_eip_sd);
  sigma_exp      ~ normal(0, 1);
  z_exp          ~ std_normal();

  // Likelihood
  for (n in 1:N) {
    real s_ctl = rho * tau_obs[n];
    real s_trt = rho * tau_obs[n];
    if (delta[n] >= 0) {
      // Cases (i) and (ii): exposure at or after infection.
      if (tau_obs[n] > delta[n]) {
        real x_arg = (tau_obs[n] - delta[n]) / s_half;
        s_trt -= B * s_half * hill_J(x_arg, nH, gl_nodes_01, gl_weights_01);
      }
    } else {
      // Case (iii): exposure before infection.
      real x_hi = (tau_obs[n] - delta[n]) / s_half;   // = (tau_obs + |delta|) / s_half
      real x_lo = -delta[n] / s_half;                 // = |delta| / s_half
      s_trt -= B * s_half * (hill_J(x_hi, nH, gl_nodes_01, gl_weights_01)
                           - hill_J(x_lo, nH, gl_nodes_01, gl_weights_01));
    }

    real logit_p_ctl = gamma_lcdf(s_ctl | shape_r, 1) - gamma_lccdf(s_ctl | shape_r, 1);
    real logit_p_trt = gamma_lcdf(s_trt | shape_r, 1) - gamma_lccdf(s_trt | shape_r, 1);

    target += binomial_logit_lpmf(pos_ctl[n] | n_ctl[n],
                                  logit_p_ctl + beta_exp[exp_id[n]]);
    target += binomial_logit_lpmf(pos_trt[n] | n_trt[n],
                                  logit_p_trt + beta_exp[exp_id[n]]);
  }
}

generated quantities {
  // Posterior predictive + per-row probabilities + pointwise log-lik.
  array[N] int pos_ctl_rep;
  array[N] int pos_trt_rep;
  vector[N] p_ctl;
  vector[N] p_trt;
  vector[N] log_lik;

  for (n in 1:N) {
    real s_ctl = rho * tau_obs[n];
    real s_trt = rho * tau_obs[n];
    if (delta[n] >= 0) {
      if (tau_obs[n] > delta[n]) {
        real x_arg = (tau_obs[n] - delta[n]) / s_half;
        s_trt -= B * s_half * hill_J(x_arg, nH, gl_nodes_01, gl_weights_01);
      }
    } else {
      real x_hi = (tau_obs[n] - delta[n]) / s_half;
      real x_lo = -delta[n] / s_half;
      s_trt -= B * s_half * (hill_J(x_hi, nH, gl_nodes_01, gl_weights_01)
                           - hill_J(x_lo, nH, gl_nodes_01, gl_weights_01));
    }

    real logit_p_ctl = gamma_lcdf(s_ctl | shape_r, 1) - gamma_lccdf(s_ctl | shape_r, 1);
    real logit_p_trt = gamma_lcdf(s_trt | shape_r, 1) - gamma_lccdf(s_trt | shape_r, 1);

    real eta_ctl = logit_p_ctl + beta_exp[exp_id[n]];
    real eta_trt = logit_p_trt + beta_exp[exp_id[n]];

    p_ctl[n] = inv_logit(eta_ctl);
    p_trt[n] = inv_logit(eta_trt);

    pos_ctl_rep[n] = binomial_rng(n_ctl[n], p_ctl[n]);
    pos_trt_rep[n] = binomial_rng(n_trt[n], p_trt[n]);

    log_lik[n] = binomial_logit_lpmf(pos_ctl[n] | n_ctl[n], eta_ctl)
               + binomial_logit_lpmf(pos_trt[n] | n_trt[n], eta_trt);
  }
}
