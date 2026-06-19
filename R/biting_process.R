#' @title Biting process
#' @description
#' This is the biting process. It results in human and mosquito infection and
#' mosquito death.
#' @param renderer the model renderer object
#' @param solvers mosquito ode solvers
#' @param models mosquito ode models
#' @param variables a list of all of the model variables
#' @param events a list of all of the model events
#' @param parameters model pararmeters
#' @param lagged_infectivity a list of LaggedValue objects with historical sums
#' of infectivity, one for every metapopulation
#' @param lagged_eir a LaggedValue class with historical EIRs
#' @param mixing_fn a function to retrieve the mixed EIR and infectivity based
#' on the other populations
#' @param mixing_index an index for this population's position in the
#' lagged_infectivity list (default: 1)
#' @param infection_outcome competing hazards object for infection rates
#' @param timestep the current timestep
#' @noRd
create_biting_process <- function(
  renderer,
  solvers,
  models,
  variables,
  events,
  parameters,
  lagged_infectivity,
  lagged_eir,
  mixing_fn = NULL,
  mixing_index = 1,
  infection_outcome
  ) {
  function(timestep) {
    # Calculate combined EIR
    age <- get_age(variables$birth$get_values(), timestep)
    bitten <- simulate_bites(
      renderer,
      solvers,
      models,
      variables,
      events,
      age,
      parameters,
      timestep,
      lagged_infectivity,
      lagged_eir,
      mixing_fn,
      mixing_index
    )
    
    simulate_infection(
      variables,
      events,
      bitten$bitten_humans,
      bitten$n_bites_per_person,
      age,
      parameters,
      timestep,
      renderer,
      infection_outcome
    )
  }
}

#' @importFrom stats rpois
simulate_bites <- function(
  renderer,
  solvers,
  models,
  variables,
  events,
  age,
  parameters,
  timestep,
  lagged_infectivity,
  lagged_eir,
  mixing_fn = NULL,
  mixing_index = 1
  ) {
  bitten_humans <- individual::Bitset$new(parameters$human_population)
  n_bites_per_person <- numeric(0)
  
  human_infectivity <- variables$infectivity$get_values()
  if (parameters$tbv) {
    human_infectivity <- account_for_tbv(
      timestep,
      human_infectivity,
      variables,
      parameters
    )
  }
  renderer$render('infectivity', mean(human_infectivity), timestep)
  
  # Calculate pi (the relative biting rate for each human)
  psi <- unique_biting_rate(age, parameters)
  zeta <- variables$zeta$get_values()
  .pi <- human_pi(zeta, psi)
  
  # Get some indices for later
  if (parameters$individual_mosquitoes) {
    infectious_index <- variables$mosquito_state$get_index_of('Im')
    susceptible_index <- variables$mosquito_state$get_index_of('Sm')
    adult_index <- variables$mosquito_state$get_index_of('NonExistent')$not(TRUE)
  }
  
  EIR <- 0
  n_bites_per_person <- rep(0, length(psi))
  
  for (s_i in seq_along(parameters$species)) {
    species_name <- parameters$species[[s_i]]
    solver_states <- solvers[[s_i]]$get_states()
    p_bitten <- prob_bitten(timestep, variables, s_i, parameters)
    Q0 <- parameters$Q0[[s_i]]
    W <- average_p_successful(p_bitten$prob_bitten_survives, .pi, Q0)
    Z <- average_p_repelled(p_bitten$prob_repelled, .pi, Q0)
    f <- blood_meal_rate(s_i, Z, parameters)
    a <- .human_blood_meal_rate(f, s_i, W, parameters)
    lambda <- effective_biting_rates(a, .pi, p_bitten)

    if (parameters$individual_mosquitoes) {
      species_index <- variables$species$get_index_of(
        parameters$species[[s_i]]
      )$and(adult_index)
      n_infectious <- calculate_infectious_individual(
        s_i,
        variables,
        infectious_index,
        adult_index,
        species_index,
        parameters
      )
    } else {
      n_infectious <- calculate_infectious_compartmental(solver_states, parameters)
    }
    
    # store the current population's EIR for later
    lagged_eir[[s_i]]$save(
      n_infectious * a,
      timestep
    )

    # lagged EIR
    if (is.null(mixing_fn)) {
      species_eir <- lagged_eir[[s_i]]$get(timestep - parameters$de)
    } else {
      species_eir <- mixing_fn(timestep=timestep)$eir[[mixing_index, s_i]]
    }

    renderer$render(paste0('EIR_', species_name), species_eir, timestep)
    EIR <- EIR + species_eir

    expected_bites <- species_eir * mean(psi)
    if (expected_bites > 0) {
      n_bites <- rpois(1, expected_bites)
      if (n_bites > 0) {
        bitten <- fast_weighted_sample(n_bites, lambda)
        bitten_humans$insert(bitten)
        renderer$render('n_bitten', bitten_humans$size(), timestep)
        if(parameters$parasite == "vivax"){
          # p.v must pass through the number of bites per person
          n_bites_per_person <- n_bites_per_person + tabulate(bitten, nbins = length(lambda))
        }
      }
    }

    lagged_infectivity$save(sum(human_infectivity * .pi), timestep)

    if (is.null(mixing_fn)) {
      infectivity <- lagged_infectivity$get(timestep - parameters$delay_gam)
    } else {
      infectivity <- mixing_fn(timestep=timestep)$inf[[mixing_index]]
    }

    foim <- calculate_foim(a, infectivity)
    renderer$render(paste0('FOIM_', species_name), foim, timestep)
    mu <- death_rate(f, W, Z, s_i, parameters)
    renderer$render(paste0('mu_', species_name), mu, timestep)
    
    if (parameters$individual_mosquitoes) {
      # update the ODE with stats for ovoposition calculations
      aquatic_mosquito_model_update(
        models[[s_i]]$.model,
        species_index$size(),
        f,
        mu
      )
      
      # update the individual mosquitoes
      susceptible_species_index <- susceptible_index$copy()$and(species_index)
      
      biting_effects_individual(
        variables,
        foim,
        events,
        s_i,
        susceptible_species_index,
        species_index,
        mu,
        parameters,
        timestep
      )
    } else {
      kernels <- compute_atn_kernels(timestep, parameters, foim, s_i)
      adult_mosquito_model_update(
        models[[s_i]]$.model,
        mu,
        foim,
        a * kernels$delta_atn * kernels$contact_factor,  # av_da = contact rate * exposure probability
        kernels$delta_atn,
        kernels$dn_atn,
        kernels$Lambda0_t,
        kernels$Lambda_i,
        kernels$rho_i,
        kernels$B_post,
        f
      )
    }
  }

  list(bitten_humans = bitten_humans, n_bites_per_person = n_bites_per_person)
}


# =================
# Utility functions
# =================

calculate_eir <- function(species, solvers, variables, parameters, timestep) {
  a <- human_blood_meal_rate(species, variables, parameters, timestep)
  infectious <- calculate_infectious(species, solvers, variables, parameters)
  infectious * a
}

effective_biting_rates <- function(a, .pi, p_bitten) {
  a * .pi * p_bitten$prob_bitten / sum(.pi * p_bitten$prob_bitten_survives)
}

calculate_infectious <- function(species, solvers, variables, parameters) {
  if (parameters$individual_mosquitoes) {
    adult_index <- variables$mosquito_state$get_index_of('NonExistent')$not(TRUE)
    return(
      calculate_infectious_individual(
        species,
        variables,
        variables$mosquito_state$get_index_of('Im'),
        adult_index,
        variables$species$get_index_of(
          parameters$species[[species]]
        )$and(adult_index),
        parameters
      )
    )
  }
  calculate_infectious_compartmental(solvers[[species]]$get_states(), parameters)
}

calculate_infectious_individual <- function(
  species,
  variables,
  infectious_index,
  adult_index,
  species_index,
  parameters
  ) {
  infectious_index$copy()$and(species_index)$size()
}

calculate_infectious_compartmental <- function(solver_states, parameters) {
  iv_idx <- iv_block_indices(parameters$deltaq, parameters$spor_len)
  max(sum(solver_states[iv_idx]), 0)
}

# Compute per-timestep ATN kernel scalars and vectors (v3 lines 773-864 + 430-467).
# Returns: delta_atn, dn_atn, Lambda0_t, Lambda_i (len deltaqp1),
#          rho_i (len deltaqp1), B_post (len spor_len).
# Lambda_i[1] is always foim (baseline; no Hill decay for unexposed compartment).
compute_atn_kernels <- function(timestep, parameters, foim, species) {
  deltaqp1 <- parameters$deltaq + 1L
  spor_len <- parameters$spor_len
  n        <- parameters$n_atn
  t0       <- parameters$t0_atn
  Q0       <- parameters$Q0_atn
  rho      <- spor_len / parameters$dem

  # --- coverage: per-event with random proportional replacement ---
  # lambda_atn is NULL when no set_bednets call has been made (no nets at all);
  # fall back to 0 (no waning) in that case.
  lambda_atn <- if (is.null(parameters$lambda_atn)) 0 else parameters$lambda_atn
  repl_factor <- vapply(seq_len(n), function(i) {
    later <- which(t0 > t0[i] & t0 <= timestep)
    prod(1 - Q0[later])
  }, numeric(1))
  Q_each <- ifelse(timestep < t0, 0,
              Q0 * exp(-lambda_atn * (timestep - t0)) * repl_factor)
  Q_t <- sum(Q_each)

  # --- per-event drug-effect decay ---
  Lambda   <- foim
  Lambda00 <- Lambda * (1 - parameters$B_max_post)   # pre- & post-infection blocking share b_max
  rho00    <- parameters$rho_frac * rho
  age      <- pmax(timestep - t0, 0)
  Lambda0_each <- ifelse(timestep < t0, Lambda,
                    Lambda - (Lambda - Lambda00) * exp(-parameters$gamma_atn * age))
  rho0_each    <- ifelse(timestep < t0, rho,
                    rho    - (rho    - rho00)    * exp(-parameters$gamma_atn * age))
  dn_each      <- ifelse(timestep < t0, 0,
                    parameters$dn0_atn * exp(-parameters$gamma_atn * age))

  # --- coverage-weighted averages ---
  Lambda0_t <- if (Q_t > 0) sum(Q_each * Lambda0_each) / Q_t else Lambda
  rho0_t    <- if (Q_t > 0) sum(Q_each * rho0_each)    / Q_t else rho
  dn_atn    <- if (Q_t > 0) sum(Q_each * dn_each)      / Q_t else 0

  # --- contact_factor: barrier-repelled mosquitoes (prob rnm) physically touch the
  # net and pick up the drug; only chemical excito-repellency (rn - rnm) prevents
  # contact entirely. So the ATN-exposure rate scales by
  #   contact_factor = (sn + rnm) / sn = (1 - rn_chem - dn) / (1 - rn - dn)
  # where sn = 1 - rn - dn (feed-and-survive probability).
  # Non-insecticidal ATN (rn0 = rnm, dn0 = 0): contact_factor = 1 / (1 - rnm).
  # Pyr-ATN: starts at (1 - rn0 - dn0 + rnm)/(1 - rn0 - dn0) for a fresh net,
  #          decays toward 1/(1 - rnm) as insecticide wanes.
  # Source rn/rnm/dn0/gamman by matching each t0_atn to its bednet schedule row.
  # Falls back to 1 (no adjustment) when no bednet schedule or no row match.
  contact_factor <- if (is.null(parameters$bednet_timesteps)) {
    1
  } else {
    bed_idx  <- match(t0, parameters$bednet_timesteps)
    idx_safe <- ifelse(!is.na(bed_idx), bed_idx, 1L)  # safe subscript; unmatched overridden below

    rn0_e  <- parameters$bednet_rn [idx_safe, species]
    rnm_e  <- parameters$bednet_rnm[idx_safe, species]
    dn0_e  <- parameters$bednet_dn0[idx_safe, species]
    gam_e  <- parameters$bednet_gamman[idx_safe]

    decay_e <- exp(-age / gam_e)                            # bednet_decay() inline
    rn_e    <- (rn0_e - rnm_e) * decay_e + rnm_e           # rn(dt): prob_repelled_bednets
    dn_e    <- dn0_e * decay_e                              # dn(dt): prob_survives_bednets
    sn_e    <- pmax(1 - rn_e - dn_e, 1e-6)                 # floor avoids divide-by-zero
    cf_each <- ifelse(!is.na(bed_idx), (sn_e + rnm_e) / sn_e, 1)

    if (Q_t > 0) sum(Q_each * cf_each) / Q_t else 1
  }

  delta_atn <- parameters$p_atn * parameters$phi_bednets[[species]] * Q_t

  # --- Bompard TRA -> field TBA transform ---
  bompard <- function(b_lab) {
    m <- parameters$m_bompard; k <- parameters$k_bompard
    a_b <- (k / (k + m))^k
    b_b <- (k / (k + m * (1 - b_lab)))^k
    (b_b - a_b) / (1 - a_b)
  }

  # --- Lambda_i: per-compartment FOI (Hill decay over ATN-exposure index) ---
  # s_days(i) = (i - 1.5) * (atn_window / deltaq): midpoint time-since-exposure in days.
  # With defaults (atn_window=10, deltaq=10) this equals (i - 1.5), unchanged.
  # Compartment 1 = unexposed baseline; s[1] = 0 guard avoids (-0.5)^nH NaN
  # (value unused — Lambda_i[1] and rho_i[1] are overwritten to baseline below).
  s    <- (seq_len(deltaqp1) - 1.5) * (parameters$atn_window / parameters$deltaq)
  s[1] <- 0
  b_lab_pre   <- (1 - Lambda0_t / Lambda) *
    (parameters$s_half_pre^parameters$nH_pre /
     (parameters$s_half_pre^parameters$nH_pre + s^parameters$nH_pre))
  b_field_pre <- if (parameters$use_bompard) bompard(b_lab_pre) else b_lab_pre
  Lambda_i    <- Lambda * (1 - b_field_pre)
  Lambda_i[1] <- Lambda   # baseline compartment always carries raw Lambda

  # --- rho_i: per-compartment EIP rate ---
  rho_i <- if (parameters$use_eip_hill) {
    rho - (rho - rho0_t) *
      (parameters$s_half_eip^parameters$nH_eip /
       (parameters$s_half_eip^parameters$nH_eip + s^parameters$nH_eip))
  } else {
    rho - (rho - rho0_t) * exp(-parameters$zeta * s)
  }
  rho_i[1] <- rho   # baseline compartment uses scalar rho (C++ ignores rho_i[1])

  # --- B_post: post-infection blocking probability ---
  t_post     <- (seq_len(spor_len) - 0.5) * parameters$dem / spor_len
  b_lab_post <- parameters$B_max_post *
    (parameters$s_half_post^parameters$nH_post /
     (parameters$s_half_post^parameters$nH_post + t_post^parameters$nH_post))
  B_post <- if (parameters$use_bompard) bompard(b_lab_post) else b_lab_post

  list(
    delta_atn      = delta_atn,
    contact_factor = contact_factor,
    dn_atn         = dn_atn,
    Lambda0_t      = Lambda0_t,
    Lambda_i       = Lambda_i,
    rho_i          = rho_i,
    B_post         = B_post
  )
}

intervention_coefficient <- function(p_bitten) {
  p_bitten$prob_bitten / sum(p_bitten$prob_bitten_survives)
}

human_pi <- function(zeta, psi) {
  (zeta * psi) / sum(zeta * psi)
}

blood_meal_rate <- function(v, z, parameters) {
  gonotrophic_cycle <- get_gonotrophic_cycle(v, parameters)
  interrupted_foraging_time <- parameters$foraging_time[[v]] / (1 - z)
  1 / (interrupted_foraging_time + gonotrophic_cycle)
}

human_blood_meal_rate <- function(species, variables, parameters, timestep) {
  age <- get_age(variables$birth$get_values(), timestep)
  psi <- unique_biting_rate(age, parameters)
  zeta <- variables$zeta$get_values()
  p_bitten <- prob_bitten(timestep, variables, species, parameters)
  .pi <- human_pi(zeta, psi)
  Q0 <- parameters$Q0[[species]]
  W <- average_p_successful(p_bitten$prob_bitten_survives, .pi, Q0)
  Z <- average_p_repelled(p_bitten$prob_repelled, .pi, Q0)
  f <- blood_meal_rate(species, Z, parameters)
  .human_blood_meal_rate(f, species, W, parameters)
}

.human_blood_meal_rate <- function(f, v, W, parameters) {
  Q <- 1 - (1 - parameters$Q0[[v]]) / W
  Q * f
}

average_p_repelled <- function(p_repelled, .pi, Q0) {
  Q0 * sum(.pi * p_repelled)
}

average_p_successful <- function(prob_bitten_survives, .pi, Q0) {
  (1 - Q0) + Q0 * sum(.pi *  prob_bitten_survives)
}

# Unique biting rate (psi) for a human of a given age
unique_biting_rate <- function(age, parameters) {
  1 - parameters$rho * exp(- age / parameters$a0)
}

#' @title Calculate the force of infection towards mosquitoes
#'
#' @param a human blood meal rate
#' @param infectivity_sum the sum of each individual's infectivity 
#' @noRd
calculate_foim <- function(a, infectivity_sum) {
  a * infectivity_sum
}
