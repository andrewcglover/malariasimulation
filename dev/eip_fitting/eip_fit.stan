// eip_fit.stan
// ATN EIP-slowdown model fit.
//
// Mechanism (for fresh nets, tau = 0 in the PDF):
//   After exposure at time Delta relative to infection, progression rate
//   through the EIP recovers from rho_0 back to rho with rate zeta:
//     rho(s) = rho - (rho - rho_0) * exp(-zeta * s),   s = time since exposure
//   so the speed deficit (rho - rho(s)) decays exponentially at rate zeta.
//
// Time-change argument:
//   EIP completion time has the same distribution as the hitting time of
//   Gamma(shape = Delta_r, rate = 1) by the integrated rate s(tau_obs; Delta).
//   Hence
//     P(EIP completed by tau_obs | Delta) = F_Gamma(s(tau_obs; Delta); Delta_r, 1).
//
//   Closed forms (with A = (rho - rho_0) / zeta):
//
//   Case (i)  Delta >= 0, tau_obs > Delta (exposure has occurred by readout):
//     s(tau_obs; Delta) = rho * tau_obs - A * (1 - exp(-zeta * (tau_obs - Delta)))
//
//   Case (ii) Delta >= 0, tau_obs <= Delta (no exposure yet by readout):
//     s(tau_obs; Delta) = rho * tau_obs
//
//   Case (iii) Delta < 0 (exposure before infection):
//     s(tau_obs; Delta) = rho * tau_obs
//                       - A * exp(zeta * Delta) * (1 - exp(-zeta * tau_obs))
//     Derivation: with exposure at Delta < 0 and infection at tau = 0, the
//     suppression has been recovering for |Delta| days before infection, so
//       integral_{0}^{tau_obs} sigma(tau - Delta) dtau
//         = (1/zeta) * (exp(-zeta*|Delta|) - exp(-zeta*(tau_obs+|Delta|)))
//         = (1/zeta) * exp(zeta*Delta) * (1 - exp(-zeta * tau_obs)),
//     and s = rho*tau_obs minus B times that integral, with B = rho - rho_0
//     and A = B/zeta. The factor exp(zeta*Delta) (with Delta < 0, so in
//     (0,1)) damps the suppression integral by however much recovery has
//     already occurred before the parasite was introduced.
//
// Likelihood: paired Binomial per row with a shared logit-scale OLRE offset
// beta_exp[exp_id] applied identically to control and 200mg arms, so it
// cancels in the within-row contrast (protecting rho_0 / zeta identifiability).

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

  // EIP prior hyperparameters (sourced from priors_atn_main.R).  Naming
  // convention: "excess" = the part of the baseline (control) EIP above
  // the biological floor; "extension" is reserved for the ATN-induced
  // lengthening of the EIP in the exposed cohort, which is governed
  // jointly by r0_frac (depressed-rate fraction) and log_zeta (recovery
  // rate back to baseline).  Reparameterise rho as
  //   rho = Delta_r / (eip_floor + eip_excess),   eip_excess = exp(log_eip_excess),
  // so the baseline mean EIP has a hard biological lower bound at eip_floor.
  real<lower=0> eip_floor;
  real log_eip_excess_mean;
  real<lower=0> log_eip_excess_sd;
  real<lower=0> r0_frac_alpha;
  real<lower=0> r0_frac_beta;
  real log_zeta_mean;
  real<lower=0> log_zeta_sd;
}

transformed data {
  real shape_r = Delta_r;
}

parameters {
  real log_eip_excess;               // log of baseline-EIP excess above floor (days)
  real<lower=0, upper=1> r0_frac;    // rho_0 / rho in (0, 1); ATN-induced
                                     // fractional slowdown at exposure
  real log_zeta;                     // log of recovery rate back to baseline (1/day)
  real<lower=0> sigma_exp;           // SD of OLRE on logit scale
  vector[K] z_exp;                   // non-centred OLRE draws
}

transformed parameters {
  real<lower=0> eip_excess   = exp(log_eip_excess);       // baseline EIP - eip_floor, days
  real<lower=0> eip_baseline = eip_floor + eip_excess;    // mean baseline (control) EIP, days
  real<lower=0> rho          = shape_r / eip_baseline;    // baseline rate (compartments/day)
  real<lower=0> zeta         = exp(log_zeta);
  real<lower=0> rho0         = r0_frac * rho;
  real A                     = (rho - rho0) / zeta;
  vector[K] beta_exp         = sigma_exp * z_exp;
}

model {
  // Priors (all hyperparameters from priors_atn_main.R via the data block).
  log_eip_excess ~ normal(log_eip_excess_mean, log_eip_excess_sd);
  r0_frac        ~ beta(r0_frac_alpha, r0_frac_beta);
  log_zeta       ~ normal(log_zeta_mean, log_zeta_sd);
  sigma_exp      ~ normal(0, 1);                          // half-normal via <lower=0>
  z_exp          ~ std_normal();

  // Likelihood
  for (n in 1:N) {
    real s_ctl = rho * tau_obs[n];
    real s_trt = rho * tau_obs[n];
    if (delta[n] >= 0) {
      // Cases (i) and (ii): exposure at or after infection.
      if (tau_obs[n] > delta[n]) {
        s_trt -= A * (1 - exp(-zeta * (tau_obs[n] - delta[n])));
      }
      // else case (ii): exposure not yet by readout, no suppression term.
    } else {
      // Case (iii): exposure before infection.
      s_trt -= A * exp(zeta * delta[n]) * (1 - exp(-zeta * tau_obs[n]));
    }

    // Logit of completion probability, numerically stable on both tails.
    real logit_p_ctl = gamma_lcdf(s_ctl | shape_r, 1) - gamma_lccdf(s_ctl | shape_r, 1);
    real logit_p_trt = gamma_lcdf(s_trt | shape_r, 1) - gamma_lccdf(s_trt | shape_r, 1);

    target += binomial_logit_lpmf(pos_ctl[n] | n_ctl[n],
                                  logit_p_ctl + beta_exp[exp_id[n]]);
    target += binomial_logit_lpmf(pos_trt[n] | n_trt[n],
                                  logit_p_trt + beta_exp[exp_id[n]]);
  }
}

generated quantities {
  // Interpretable derived quantities (eip_baseline already a transformed
  // parameter; expose the recovery half-life here for posterior summaries).
  real zeta_halflife = log(2) / zeta;  // recovery half-life in days

  // Posterior predictive + completion probabilities + pointwise log-lik.
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
        s_trt -= A * (1 - exp(-zeta * (tau_obs[n] - delta[n])));
      }
    } else {
      s_trt -= A * exp(zeta * delta[n]) * (1 - exp(-zeta * tau_obs[n]));
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
