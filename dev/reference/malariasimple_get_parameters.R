#' @title Create parameters for malariasimple
#' @description Helper function to provide defaults for most necessary values required for running the malariasimple model
#'
#' @param stochastic Boolean variable. Set to false for deterministic simulation
#' @param n_days Number of days for which the simulation will run.
#' @param human_pop Size of human population (count)
#' @param tsd Number of time-steps per day. Fewer is faster, more better approximates the continuous solution.
#' @param age_vector Lower bound of each age category. See function for default.
#' @param biting_groups Number of biting heterogeneity groups.
#' @param lag_rates Number of sub-compartments within FOI and FOIv which approximate delay-differential equation. Higher values are a closer approximation, but computationally more expensive.

#' @param eta Death rate for exponential population distribution, i.e. 1/Mean Population Age
#' @param rho Age-dependent biting parameter
#' @param a0 Age-dependent biting parameter
#' @param sigma2 Variance of the log heterogeneity in biting rates
#' @param max_age Upper boundary of oldest age category
#'
#' @param rA Rate of leaving asymptomatic infection
#' @param rT Rate of leaving treatment
#' @param rD Rate of leaving clinical disease
#' @param rU Rate of recovering from subpatent infection
#' @param rP Rate of leaving prophylaxis
#'
#' @param dE Latent period of human infection
#' @param delayGam Lag from parasites to infectious gametocytes
#'
#' @param cD Untreated disease contribution to infectiousness
#' @param cT Treated disease contribution to infectiousness
#' @param cU Subpatent disease contribution to infectiousness
#' @param gamma1 Parameter for infectiousness of state A
#' @param d1 Minimum probability due to maximum immunity
#' @param dID Inverse of decay rate
#' @param ID0 Scale parameter
#' @param kD Shape parameter
#' @param uD Duration in which immunity is not boosted
#' @param aD Scale parameter relating age to immunity
#' @param fD0 Time-scale at which immunity changes with age
#' @param gammaD Shape parameter relating age to immunity
#' @param alphaA PCR detection probability parameters state A
#' @param alphaU PCR detection probability parameters state U
#' @param b0 Maximum probability due to no immunity
#' @param b1 Maximum relative reduction due to immunity
#' @param dB Inverse of decay rate
#' @param IB0 Scale parameter
#' @param kB Shape parameter
#' @param uB Duration in which immunity is not boosted
#' @param phi0 Maximum probability due to no immunity
#' @param phi1 Maximum relative reduction due to immunity
#' @param dCA Inverse of decay rate
#' @param IC0 Scale parameter
#' @param kC Shape parameter
#' @param uCA Duration in which immunity is not boosted
#' @param PM New-born immunity relative to mothers
#' @param dCM Inverse of decay rate of maternal immunity
#' @param delayMos Extrinsic incubation period
#' @param foraging_time Duration of host seeking, assumed to be constant between species
#' @param gonotrophic_cycle Duration of mosquito resting after feed
#' @param mum Daily mortality of adult mosquitoes
#' @param Q0 Anthrophagy probability
#' @param chi Endophily probability
#' @param phi_bednets Percentage of bites indoors and in bed
#' @param phi_indoors Percentage of bites indoors
#' @param muEL Per capita daily mortality rate of early stage larvae (low density)
#' @param muLL Per capita daily mortality rate of late stage larvae (low density)
#' @param muPL Per capita daily mortality rate of pupae
#' @param dEL Development time of early stage larvae
#' @param dLL Development time of late stage larvae
#' @param dPL Development time of pupae
#' @param gammaL Relative effect of density dependence on late instars relative to early instars
#' @param betaL Number of eggs laid per day per mosquito
#' @param daily_ft Daily vector of percentage of population that gets treated. A scalar value is permitted and assumes constant ft.
#' @param clin_inc_rendering_min_ages Vector of values (or singe value) of lower age boundaries for clinical incidence output (days)
#' @param clin_inc_rendering_max_ages Vector of values (or singe value) of upper age boundaries for clinical incidence output (days)
#' @param prevalence_rendering_min_ages Vector of values (or singe value) of lower age boundaries for prevalence output (days)
#' @param prevalence_rendering_max_ages Vector of values (or singe value) of upper age boundaries for prevalence output (days)

#' @param ... Additional arguments
#' @examples
#' #Output default model for clinical incidence in range [0-10] and 2+
#' params <- get_parameters(n_days = 500,
#'                          clin_inc_rendering_min_ages = c(0,2*365),
#'                          clin_inc_rendering_max_ages = c(10*365,Inf)) |>
#'           set_equilibrium(init_EIR = 10)
#' @export
get_parameters <- function(
    ##Accuracy/speed trade off parameters
    stochastic = FALSE,
    n_days = 100,
    human_pop = 100000,
    tsd = 4, #Time steps per day
    age_vector =  c(0,0.25,0.5,1,1.5,2,2.5,3,3.5,4,5,6,7,8.5,10,20,30,40,60,80)*365, #Default
    biting_groups = 3,
    lag_rates = 10, #Number of sub-compartments within FOI and FOIv which approximate delay-differential equation. Higher values are a closer approximation, but computationally more expensive
    # age, heterogeneity in exposure,
    eta = 1/(21*365),
    rho = 0.85,
    a0 = 2920,
    sigma2 = 1.67,
    max_age = 100*365,
    #  rate of leaving infection states
    rA = 1/195,
    rT = 0.2,
    rD = 0.2,
    rU = 1/110.299,
    rP = 1/15,
    #  human latent period and time lag from asexual parasites to
    dE  = 12,
    delayGam = 12.5,
    # human infectiousness to mosquitoes
    cD  = 0.0676909,
    cT  =  0.322 * cD,
    cU  = 0.006203,
    gamma1  = 1.82425,
    #  Immunity reducing probability of detection
    d1 = 0.160527,
    dID = 3650,
    ID0 = 1.577533,
    kD = 0.476614,
    uD = 9.44512,
    aD = 8001.99,
    fD0 = 0.007055,
    gammaD = 4.8183,
    alphaA = 0.75735,
    alphaU = 0.185624,
    # Immunity reducing probability of infection
    b0 = 0.590076,
    b1 = 0.5,
    dB = 3650,
    IB0 = 43.8787,
    kB = 2.15506,
    uB = 7.19919,
    # Immunity reducing probability of clinical disease
    phi0 = 0.791666,
    phi1 = 0.000737,
    dCA = 10950,
    IC0 = 18.02366,
    kC = 2.36949,
    uCA = 6.06349,
    PM = 0.774368,
    dCM = 67.6952,
    # entomological parameters
    delayMos = 10,
    foraging_time = 0.69,
    gonotrophic_cycle = 2.31,
    mum = 0.132,
    Q0 = 0.92,
    chi = 0.86,
    phi_bednets = 0.85,
    phi_indoors = 0.9,
    # Erlang EIP
    spor_len = 10,
    # larval parameters daily density dependent mortality rate of egg
    muEL = 0.0338,
    muLL = 0.0348,
    muPL = 0.249,
    dEL = 6.64,
    dLL = 3.72,
    dPL = 0.643,
    gammaL = 13.25,
    betaL = 21.2,
    # ATN parameters
    t0_atn     = 0,  # distribution time (days)
    Q0_atn     = 0,  # initial coverage right after distribution
    lambda_atn = 0.0017, # mean loss rate (1 / mean duration)
    p_atn      = 0.9,  # prob(exposure | hits a net)
    deltaq     = 10,   # number of ATN-exposed compartments
    gamma_atn  = 0.000719, #log(2) / (2.64 * 365)
    xi         = 0.3,
    zeta       = 0.433217,#0.1732868,#0.1386294,
    Lambda00sf = 0.0019  ,
    #rho00      = 0.12, # superseeded by rho_frac
    dn0_atn    = 0,
    # ATN pre-infection blocking Hill parameters
    s_half_pre = 1.1539  ,
    nH_pre     = 4.0229  ,
    # ATN post-infection blocking Hill parameters
    B_max_post  = 0.9981  ,
    s_half_post = 0.6819  ,
    nH_post     = 14.1709 ,
    # ATN EIP lengthening Hill parameters
    s_half_eip = 2.6735  ,
    nH_eip     = 4.8142 ,
    rho_frac    = 0.1340  ,
    # ATN Bompard parameters
    m_bompard = 0.000157,
    k_bompard = 4.95e-6,
    # ATN logic
    use_bompard = TRUE,
    use_eip_hill = TRUE,

    # intervention parameters
    daily_ft = 0,
    clin_inc_rendering_min_ages = NULL,
    clin_inc_rendering_max_ages = NULL,
    prevalence_rendering_min_ages = NULL, #Default = 2*365
    prevalence_rendering_max_ages = NULL, #Default = 10*365
    ...

){
  # set up param list
  params <- list()

  # catch extra params and place in list
  extra_param_list <- list(...)
  if(length(extra_param_list)>0){
    if(is.list(extra_param_list[[1]])){
      extra_param_list <- extra_param_list[[1]]
    }
  }

  #Sanity Checks
  if(age_vector[1] != 0){stop(message("age_vector must start from zero"))}
  if(length(age_vector) < 2){stop(message("age_vector must have at least two categories"))}
  if(max(age_vector) <= 20*365){stop(message("At least one age category must be over 20 years (7300 days)"))}
  if(!is.numeric(n_days)){stop(message("n_days must be a numeric value"))}
  if(!is.numeric(tsd)){stop(message("tsd must be a numeric value"))}
  if(!is.numeric(biting_groups)){stop(message("biting_groups must be a numeric value"))}
  if(!is.numeric(lag_rates)){stop(message("biting_groups must be a numeric value"))}
  if(n_days %% 1 != 0 | n_days < 1){stop(message("n_days must be a positive integer value"))}
  if(tsd %% 1 != 0 | tsd < 1){stop(message("tsd must be a positive integer value"))}
  if(biting_groups %% 1 != 0 | biting_groups < 1){stop(message("biting_groups must be a positive integer value"))}
  if(lag_rates %% 1 != 0| lag_rates < 1){stop(message("lag_rates must be a positive integer value"))}
  if(!is.numeric(daily_ft)){stop(message("daily_ft may only contain numeric values between 0 and 1"))}
  if(any(daily_ft < 0 | daily_ft > 1)){stop(message("daily_ft may only contain numeric values between 0 and 1"))}
  if(length(daily_ft) < n_days && length(daily_ft) != 1){stop(message("daily_ft must be either scaler or at least as long as n_days"))}

  # ATN distribution event checks. t0_atn and Q0_atn may be passed as
  # scalars (a single distribution event) or as equal-length vectors
  # (multiple events); each (t0_atn, Q0_atn) pair defines one ATN event.
  if(!is.numeric(t0_atn) | !is.numeric(Q0_atn)){stop(message("t0_atn and Q0_atn must be numeric"))}
  if(length(t0_atn) != length(Q0_atn)){stop(message("t0_atn and Q0_atn must be equal in length"))}
  if(any(Q0_atn < 0 | Q0_atn > 1)){stop(message("Q0_atn may only contain values between 0 and 1"))}
  if(any(t0_atn %% 1 != 0)){stop(message("t0_atn must contain integer values"))}
  if(any(t0_atn < 0)){stop(message("t0_atn values must be non-negative"))}
  if(any(diff(t0_atn) < 1)){stop(message("t0_atn values must be unique and in chronological order"))}
  if(any(t0_atn > n_days)){stop(message("t0_atn values cannot exceed n_days"))}

  ###########################################
  # Define parameters
  ###########################################
  ## DEFAULT PARAMS
  params$n_days <- n_days
  params$human_pop <- human_pop
  params$stochastic <- stochastic
  params$tsd <- tsd
  params$n_ts <- n_days*tsd

  # duration of year
  params$DY <- 365

  # age, heterogeneity in exposure
  params$age_vector <- age_vector
  params$biting_groups <- biting_groups
  params$eta <- eta
  params$rho <- rho
  params$a0 <- a0
  params$sigma2 <- sigma2
  params$max_age <- max_age
  na <- as.integer(length(age_vector))  # number of age groups
  nh <- as.integer(biting_groups)
  params$na <- na
  params$nh <- nh

  if(length(daily_ft) == 1) daily_ft <- rep(daily_ft, (n_days + 1))
  daily_ft <- daily_ft[1:n_days]
  params$daily_ft <- c(daily_ft[1], daily_ft)

  # rate of leaving infection states
  params$rA <- rA
  params$rT <- rT
  params$rD <- rD
  params$rU <- rU
  params$rP <- rP

  # human latent period and time lag from asexual parasites to
  # infectiousness
  params$dE <- dE
  params$delayGam <- delayGam

  # infectiousness to mosquitoes
  params$cD <- cD
  params$cT <- cT
  params$cU <- cU
  params$gamma1 <- gamma1

  # Immunity reducing probability of detection
  params$d1 <- d1
  params$dID <- dID
  params$ID0 <- ID0
  params$kD <- kD
  params$uD <- uD
  params$aD <- aD
  params$fD0 <- fD0
  params$gammaD <- gammaD

  # PCR prevalence parameters
  params$alphaA <- alphaA
  params$alphaU <- alphaU

  # anti-infection immunity
  params$b0 <- b0
  params$b1 <- b1
  params$dB <- dB
  params$IB0 <- IB0
  params$kB <- kB
  params$uB <- uB

  # clinical immunity
  params$phi0 <- phi0
  params$phi1 <- phi1
  params$dCA <- dCA
  params$IC0 <- IC0
  params$kC <- kC
  params$uCA <- uCA
  params$PM <- PM
  params$dCM <- dCM

  # entomological parameters
  params$delayMos <- delayMos
  params$foraging_time <- foraging_time
  params$gonotrophic_cycle <- gonotrophic_cycle
  params$mum <- mum
  params$Q0 <- Q0
  params$phi_bednets <- phi_bednets
  params$phi_indoors <- phi_indoors
  params$fv0 <- 1 / (foraging_time + gonotrophic_cycle)
  params$av0 <- Q0 * params$fv0 # daily feeeding rate on humans
  params$Surv0 <- exp(-mum * delayMos) # probability of surviving incubation period
  params$p10 <- exp(-mum * foraging_time)  # probability of surviving one feeding cycle
  params$p2 <- exp(-mum * gonotrophic_cycle)  # probability of surviving one resting cycle

  # Erlang EIP
  params$spor_len <- spor_len # number of latently-infected mosquito compartments

  # ATN parameters
  params$t0_atn     <- t0_atn  # distribution time (days)
  params$Q0_atn     <- Q0_atn  # initial coverage right after distribution
  params$n_atn      <- length(t0_atn) # number of ATN distribution events (= length of t0_atn / Q0_atn)
  params$lambda_atn <- lambda_atn # mean loss rate (1 / mean duration)
  params$p_atn      <- p_atn  # prob(exposure | hits a net)
  params$deltaq     <- deltaq # number of ATN-exposed compartments
  params$deltaqp1   <- deltaq + 1
  params$gamma_atn  <- gamma_atn
  params$xi         <- xi
  params$zeta       <- zeta
  params$Lambda00sf <- Lambda00sf
  #params$rho00     <- rho00
  params$dn0_atn    <- dn0_atn # initial probability of mortality from an ATN

  params$s_half_pre <- s_half_pre
  params$nH_pre     <- nH_pre
  # ATN post-infection blocking Hill parameters
  params$B_max_post  <- B_max_post
  params$s_half_post <- s_half_post
  params$nH_post     <- nH_post
  # ATN EIP lengthening Hill parameters
  params$s_half_eip <- s_half_eip
  params$nH_eip     <- nH_eip
  params$rho_frac    <- rho_frac
  # ATN Bomard parameters
  params$m_bompard <- m_bompard
  params$k_bompard <- k_bompard
  # ATN logical parameters
  params$use_bompard <- use_bompard
  params$use_eip_hill <- use_eip_hill

  # larval parameters
  params$muEL <- muEL
  params$muLL <- muLL
  params$muPL <- muPL
  params$dEL <- dEL
  params$dLL <- dLL
  params$dPL <- dPL
  params$gammaL <- gammaL
  params$betaL <- betaL
  # {White et al. 2011 Parasites and Vectors}
  params$eov <- betaL/mum * (exp(mum/params$fv0) - 1)
  params$b_lambda <- (gammaL * muLL/muEL - dEL/dLL + (gammaL - 1) * muLL * dEL)
  params$lambda <- -0.5 * params$b_lambda +
    sqrt(0.25 * params$b_lambda^2 + gammaL * betaL * muLL * dEL/(2 * muEL * mum * dLL * (1 + dPL * muPL)))


  #Additional parameters for dust model
  params$lag_rates <- lag_rates
  params$lag_ratesMos <- lag_rates
  params$dt <- 1/tsd

  #Set age rendering parameters for prevalence and clinical incidence estimates
  params <- age_rendering(params,
                          prevalence_rendering_max_ages,prevalence_rendering_min_ages,
                          clin_inc_rendering_max_ages,clin_inc_rendering_min_ages)

  ###########################################
  # Extras
  ###########################################
  #Track which parameters have been set so far.
  params$smc_set <- 0
  params$itn_set <- 0
  params$equilibrium_set <- 0
  params$seasonality_set <- 0

  ##Default parameters to be potentially overriden
  params$EIP_tempsens <- 0 #odin.dust likes numerical variable types
  params$mu_tempsens <- 0

  #Check that none of the spare parameters in the extra
  if(sum(!is.na(match(names(extra_param_list),names(params))))!=0){

    stop (message(cat("Extra params in ... share names with default param names. Please check:\n",
                      names(extra_param_list)[!is.na(match(names(extra_param_list),names(params)))]
    )
    ))
  }
  return(append(params,extra_param_list))
}

max_age_to_index <- function(max_age,age_vector){
  if(is.infinite(max_age)){
    max_age_index <- length(age_vector)
  } else if(!max_age %in% age_vector) {
    stop(message("rendering_max_ages must correspond to an age category boundary (or Inf)"))
  } else {
    max_age_index <- which(age_vector == max_age) - 1
  }
  return(max_age_index)
}

age_rendering <- function(params,
                          prevalence_rendering_max_ages,
                          prevalence_rendering_min_ages,
                          clin_inc_rendering_max_ages,
                          clin_inc_rendering_min_ages){
  #Set default values
  default_clin_inc_rendering_min_ages <- 0*365
  default_clin_inc_rendering_max_ages <- Inf
  default_prevalence_rendering_min_ages <- 2*365
  default_prevalence_rendering_max_ages <- 10*365

  age_vector <- params$age_vector

  #If NULL, try to apply default values. Else cover all age ranges
  if(is.null(clin_inc_rendering_min_ages)){
    if(default_clin_inc_rendering_min_ages %in% age_vector){
      clin_inc_rendering_min_ages <- default_clin_inc_rendering_min_ages
    } else {
      clin_inc_rendering_min_ages <- 0
    }
  }

  if(is.null(prevalence_rendering_min_ages)){
    if(default_prevalence_rendering_min_ages %in% age_vector){
      prevalence_rendering_min_ages <- default_prevalence_rendering_min_ages
    } else {
      prevalence_rendering_min_ages <- 0
    }
  }

  if(is.null(clin_inc_rendering_max_ages)){
    if(default_clin_inc_rendering_max_ages %in% age_vector){
      clin_inc_rendering_max_ages <- default_clin_inc_rendering_max_ages
    } else {
      clin_inc_rendering_max_ages <- Inf
    }
  }

  if(is.null(prevalence_rendering_max_ages)){
    if(default_prevalence_rendering_max_ages %in% age_vector){
      prevalence_rendering_max_ages <- default_prevalence_rendering_max_ages
    } else {
      prevalence_rendering_max_ages <- Inf
    }
  }

  ## Check rendering ages are valid.
  if(length(prevalence_rendering_max_ages) != length(prevalence_rendering_min_ages)){
    stop(message("prevalence rendering min age and max age must be equal length"))
    }
  if(length(clin_inc_rendering_max_ages) != length(clin_inc_rendering_min_ages)){
    stop(message("clinical incidence rendering min age and max age must be equal length"))
    }
  if(sum(!prevalence_rendering_min_ages %in% age_vector) != 0){
    stop(message("prevalence_rendering_min_ages must correspond to boundaries specified in age_vector. \n
                 See params$age_vector for default values"))
  }
  if(sum(!clin_inc_rendering_min_ages %in% age_vector) != 0){
    stop(message("clin_inc_rendering_min_ages must correspond to boundaries specified in age_vector. \n
                 See params$age_vector for default values"))
  }

  # Convert minimum ages (in days) to index of age vector
  params$min_age_prev <- match(prevalence_rendering_min_ages, age_vector) |> as.integer()
  params$min_age_inc <- match(clin_inc_rendering_min_ages,age_vector) |> as.integer()

  # Convert maximum ages (in days) in index of age vector. Age vector describes lower age bracket, hence maximum ages takes preceding index
  params$max_age_prev <- sapply(prevalence_rendering_max_ages, function(x) max_age_to_index(x,age_vector)) |> as.integer()
  params$max_age_inc <- sapply(clin_inc_rendering_max_ages, function(x) max_age_to_index(x,age_vector)) |> as.integer()

  # Prepare inputs for odin.dust model
  params$prev_dim <- length(prevalence_rendering_max_ages)
  params$inc_dim <- length(clin_inc_rendering_max_ages)

  # Prepare inputs for post processing
  params$prevalence_rendering_max_ages <- prevalence_rendering_max_ages
  params$clin_inc_rendering_max_ages <- clin_inc_rendering_max_ages

  return(params)
}
