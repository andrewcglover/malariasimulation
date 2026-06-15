## DECLARE TIME STEPS
n_days <- parameter()

dim(days) <- n_days + 1
days <- parameter()

dim(daily_rain_input) <- n_days +1
daily_rain_input <- parameter()
rain_input <- interpolate(days, daily_rain_input, "linear")

dim(daily_ft) <- n_days + 1
daily_ft <- parameter()
ft <- interpolate(days, daily_ft, "linear")

## MODEL VARIABLES
##------------------------------------------------------------------------------

na <- parameter()
nh <- parameter()

##------------------------------------------------------------------------------
##################
## HUMAN STATES ##
##################
##------------------------------------------------------------------------------
# Human states as specified in full transmission model
# http://journals.plos.org/plosmedicine/article?id=10.1371%2Fjournal.pmed.1000324
# http://www.nature.com/articles/ncomms4136

#--------- DEFINE PARAMETERS FOR HUMAN TRANSITIONS ---------
eta <- parameter()
age_rate <- parameter()
dim(age_rate) <- na
het_wt <- parameter()
dim(het_wt) <- nh
rA <- parameter()
rT <- parameter()
rD <- parameter()
rU <- parameter()
rP <- parameter()

#----------- INITIALISE HUMAN COMPARTMENTS -----------------
# S - SUSCEPTIBLE
init_S <- parameter()
dim(init_S) <- c(na,nh,num_int)
initial(S[,,]) <- init_S[i,j,k]
dim(S) <- c(na,nh,num_int)

# T- SUCCESSFULLY TREATED
init_T <- parameter()
dim(init_T) <- c(na,nh,num_int)
initial(T[,,]) <- init_T[i,j,k]
dim(T) <- c(na,nh,num_int)

# D - CLEAR DISEASE
init_D <- parameter()
dim(init_D) <- c(na,nh,num_int)
initial(D[,,]) <- init_D[i,j,k]
dim(D) <- c(na,nh,num_int)

# A - ASYMPTOMATIC DISEASE
init_A <- parameter()
dim(init_A) <- c(na,nh,num_int)
initial(A[,,]) <- init_A[i,j,k]
dim(A) <- c(na,nh,num_int)

# U - SUBPATENT DISEASE
init_U <- parameter()
dim(init_U) <- c(na,nh,num_int)
initial(U[,,]) <- init_U[i,j,k]
dim(U) <- c(na,nh,num_int)

# P - PROPHYLAXIS
init_P <- parameter()
dim(init_P) <- c(na,nh,num_int)
initial(P[,,]) <- init_P[i,j,k]
dim(P) <- c(na,nh,num_int)

#-------------------- DEFINE HUMAN TRANSITION EVENTS -------------------------
#### BIRTHS ####
dim(births) <- c(1,nh,num_int)
births[1,,] <- H*dt*cov[k]*eta*het_wt[j]


#### TRANSITIONS FROM S #####
# Rates of transition
dim(ST_rate) <- c(na,nh,num_int)
dim(SD_rate) <- c(na,nh,num_int)
dim(SA_rate) <- c(na,nh,num_int)

ST_rate[,,] <- ft*phi[i,j,k]*FOI_smc[i,j,k]
SD_rate[,,] <- (1-ft)*phi[i,j,k]*FOI_smc[i,j,k]
SA_rate[,,] <- (1-phi[i,j,k])*FOI_smc[i,j,k]


# Number of humans transitioning per timestep
dim(ST_trans) <- c(na,nh,num_int)
dim(SD_trans) <- c(na,nh,num_int)
dim(SA_trans) <- c(na,nh,num_int)
dim(S_death) <- c(na,nh,num_int)
dim(S_age) <- c(na,nh,num_int)

ST_trans[,,] <- if(S[i,j,k]*dt*ST_rate[i,j,k] < 0) 0 else S[i,j,k]*dt*ST_rate[i,j,k]
SD_trans[,,] <- if(S[i,j,k]*dt*SD_rate[i,j,k] < 0) 0 else S[i,j,k]*dt*SD_rate[i,j,k]
SA_trans[,,] <- if(S[i,j,k]*dt*SA_rate[i,j,k] < 0) 0 else S[i,j,k]*dt*SA_rate[i,j,k]
S_death[,,] <- if(S[i,j,k]*dt*eta < 0) 0 else S[i,j,k]*dt*eta
S_age[,,] <- if(S[i,j,k]*dt*age_rate[i] < 0) 0 else S[i,j,k]*dt*age_rate[i]

#### TRANSITIONS FROM T #####
dim(TP_trans) <- c(na,nh,num_int)
dim(T_death) <- c(na,nh,num_int)
dim(T_age) <- c(na,nh,num_int)

TP_trans[,,] <- if(T[i,j,k]*dt*rT < 0) 0 else T[i,j,k]*dt*rT
T_death[,,] <- if(T[i,j,k]*dt*eta < 0) 0 else T[i,j,k]*dt*eta
T_age[,,] <- if(T[i,j,k]*dt*age_rate[i] < 0) 0 else T[i,j,k]*dt*age_rate[i]

#### TRANSITIONS FROM D #####
dim(DA_trans) <- c(na,nh,num_int)
dim(DS_trans) <- c(na,nh,num_int)
dim(D_death) <- c(na,nh,num_int)
dim(D_age) <- c(na,nh,num_int)

DA_trans[,,] <- if(D[i,j,k]*dt*rD < 0) 0 else D[i,j,k]*dt*rD
DS_trans[,,] <- if(D[i,j,k]*alpha_smc_array[i,j,k] < 0) 0 else D[i,j,k]*alpha_smc_array[i,j,k]
D_death[,,] <- if(D[i,j,k]*dt*eta < 0) 0 else D[i,j,k]*dt*eta
D_age[,,] <- if(D[i,j,k]*dt*age_rate[i] < 0) 0 else D[i,j,k]*dt*age_rate[i]

#### TRANSITIONS FROM A #####
#Rates
dim(AT_rate) <- c(na,nh,num_int)
dim(AD_rate) <- c(na,nh,num_int)

AT_rate[,,] <- phi[i,j,k] * FOI_smc[i,j,k] * ft
AD_rate[,,] <- phi[i,j,k] * FOI_smc[i,j,k] * (1-ft)

# Number of humans transitioning per timestep
dim(AU_trans) <- c(na,nh,num_int)
dim(AT_trans) <- c(na,nh,num_int)
dim(AD_trans) <- c(na,nh,num_int)
dim(AS_trans) <- c(na,nh,num_int) #SMC Infection Clearing
dim(A_age) <- c(na,nh,num_int)
dim(A_death) <- c(na,nh,num_int)

AU_trans[,,] <- if(A[i,j,k]*dt*rA < 0) 0 else A[i,j,k]*dt*rA
AT_trans[,,] <- if(A[i,j,k]*dt*AT_rate[i,j,k] < 0) 0 else A[i,j,k]*dt*AT_rate[i,j,k]
AD_trans[,,] <- if(A[i,j,k]*dt*AD_rate[i,j,k] < 0) 0 else A[i,j,k]*dt*AD_rate[i,j,k]
AS_trans[,,] <- if(A[i,j,k]*alpha_smc_array[i,j,k] < 0) 0 else A[i,j,k]*alpha_smc_array[i,j,k] #smc Infection Clearing
A_death[,,] <- if(A[i,j,k]*dt*eta < 0) 0 else A[i,j,k]*dt*eta
A_age[,,] <- if(A[i,j,k]*dt*age_rate[i] < 0) 0 else A[i,j,k]*dt*age_rate[i]

#### TRANSITIONS FROM U #####
#Rates
dim(UA_rate) <- c(na,nh,num_int)
dim(UD_rate) <- c(na,nh,num_int)
dim(UT_rate) <- c(na,nh,num_int)

UA_rate[,,] <- (1-phi[i,j,k]) * FOI_smc[i,j,k]
UD_rate[,,] <- phi[i,j,k] * (1-ft) * FOI_smc[i,j,k]
UT_rate[,,] <- phi[i,j,k] * ft * FOI_smc[i,j,k]

# Number of humans transitioning per timestep
dim(UA_trans) <- c(na,nh,num_int)
dim(UD_trans) <- c(na,nh,num_int)
dim(UT_trans) <- c(na,nh,num_int)
dim(US_trans) <- c(na,nh,num_int)
dim(US_trans_SMC) <- c(na,nh,num_int)
dim(U_age) <- c(na,nh,num_int)
dim(U_death) <- c(na,nh,num_int)

US_trans[,,] <- if(U[i,j,k]*dt*rU < 0) 0 else U[i,j,k]*dt*rU
US_trans_SMC[,,] <- if(U[i,j,k]*alpha_smc_array[i,j,k] < 0) 0 else U[i,j,k]*alpha_smc_array[i,j,k]
UA_trans[,,] <- if(U[i,j,k]*dt*UA_rate[i,j,k] < 0) 0 else U[i,j,k]*dt*UA_rate[i,j,k]
UD_trans[,,] <- if(U[i,j,k]*dt*UD_rate[i,j,k] < 0) 0 else U[i,j,k]*dt*UD_rate[i,j,k]
UT_trans[,,] <- if(U[i,j,k]*dt*UT_rate[i,j,k] < 0) 0 else U[i,j,k]*dt*UT_rate[i,j,k]
U_age[,,] <- if(U[i,j,k]*dt*age_rate[i] < 0) 0 else U[i,j,k]*dt*age_rate[i]
U_death[,,] <- if(U[i,j,k]*dt*eta < 0) 0 else U[i,j,k]*dt*eta

#### TRANSITIONS FROM P #####
dim(PS_trans) <- c(na,nh,num_int)
dim(P_death) <- c(na,nh,num_int)
dim(P_age) <- c(na,nh,num_int)

PS_trans[,,] <-if(P[i,j,k]*dt*rP < 0) 0 else P[i,j,k]*dt*rP
P_death[,,] <- if(P[i,j,k]*dt*eta < 0) 0 else P[i,j,k]*dt*eta
P_age[,,] <- if(P[i,j,k]*dt*age_rate[i] < 0) 0 else P[i,j,k]*dt*age_rate[i]

#------------------------------ TRANSITIONS ----------------------------------
update(S[1,1:nh,1:num_int]) <- S[i,j,k] + PS_trans[i,j,k] + US_trans[i,j,k] + US_trans_SMC[i,j,k] + AS_trans[i,j,k] + DS_trans[i,j,k] - ST_trans[i,j,k] - SD_trans[i,j,k] - SA_trans[i,j,k] - S_death[i,j,k] - S_age[i,j,k] + births[1,j,k]
update(S[2:na,1:nh,1:num_int]) <- S[i,j,k] + PS_trans[i,j,k] + US_trans[i,j,k] + US_trans_SMC[i,j,k] + AS_trans[i,j,k] + DS_trans[i,j,k] - ST_trans[i,j,k] - SD_trans[i,j,k] - SA_trans[i,j,k] - S_death[i,j,k] - S_age[i,j,k] + S_age[i-1,j,k]

update(T[1,1:nh,1:num_int]) <- T[i,j,k] + AT_trans[i,j,k] + ST_trans[i,j,k] + UT_trans[i,j,k] - TP_trans[i,j,k] - T_age[i,j,k] - T_death[i,j,k]
update(T[2:na,1:nh,1:num_int]) <- T[i,j,k] + AT_trans[i,j,k] + ST_trans[i,j,k] + UT_trans[i,j,k] - TP_trans[i,j,k] - T_age[i,j,k] - T_death[i,j,k] + T_age[i-1,j,k]

update(D[1,1:nh,1:num_int]) <- D[i,j,k] + SD_trans[i,j,k] + AD_trans[i,j,k] + UD_trans[i,j,k] - DA_trans[i,j,k] - DS_trans[i,j,k] - D_death[i,j,k] - D_age[i,j,k]
update(D[2:na,1:nh,1:num_int]) <- D[i,j,k] + SD_trans[i,j,k] + AD_trans[i,j,k] + UD_trans[i,j,k] - DA_trans[i,j,k] - DS_trans[i,j,k] - D_death[i,j,k] - D_age[i,j,k] + D_age[i-1,j,k]

update(A[1,1:nh,1:num_int]) <- A[i,j,k] + SA_trans[i,j,k] + DA_trans[i,j,k] + UA_trans[i,j,k] - AT_trans[i,j,k] - AD_trans[i,j,k] - AU_trans[i,j,k] - AS_trans[i,j,k] - A_death[i,j,k] - A_age[i,j,k]
update(A[2:na,1:nh,1:num_int]) <- A[i,j,k] + SA_trans[i,j,k] + DA_trans[i,j,k] + UA_trans[i,j,k] - AT_trans[i,j,k] - AD_trans[i,j,k] - AU_trans[i,j,k] - AS_trans[i,j,k] - A_death[i,j,k] - A_age[i,j,k] + A_age[i-1,j,k]

update(U[1,1:nh,1:num_int]) <- U[i,j,k] + AU_trans[i,j,k] - UD_trans[i,j,k] - UT_trans[i,j,k] - US_trans[i,j,k] - US_trans_SMC[i,j,k] - UA_trans[i,j,k] - U_age[i,j,k] - U_death[i,j,k]
update(U[2:na,1:nh,1:num_int]) <- U[i,j,k] + AU_trans[i,j,k] - UD_trans[i,j,k] - UT_trans[i,j,k] - US_trans[i,j,k] - US_trans_SMC[i,j,k] - UA_trans[i,j,k] - U_age[i,j,k] - U_death[i,j,k] + U_age[i-1,j,k]

update(P[1,1:nh,1:num_int]) <- P[i,j,k] + TP_trans[i,j,k] - PS_trans[i,j,k] - P_death[i,j,k] - P_age[i,j,k]
update(P[2:na,1:nh,1:num_int]) <- P[i,j,k] + TP_trans[i,j,k] - PS_trans[i,j,k] - P_death[i,j,k] - P_age[i,j,k] + P_age[i-1,j,k]

##------------------------------------------------------------------------------
##########################################
## KEY INFECTION PARAMETERS (FOI & EIR) ##
##########################################
##------------------------------------------------------------------------------

#----------FORCE OF INFECTION-----------------
# Force of infection, depends on level of infection blocking immunity
dim(FOI_lag) <- c(na,nh,num_int)
FOI_lag[1:na, 1:nh, 1:num_int] <- EIR[i,j,k] * (if(IB[i,j,k]==0) b0 else b[i,j,k])

# Current FOI depends on humans that have been through the latent period
dE <- parameter()
lag_rates <- parameter(type = "integer")

dim(FOI_eq) <- c(na,nh)
FOI_eq <- parameter()

dim(FOI_XL) <- c(na,nh,num_int,lag_rates)
initial(FOI_XL[,,,]) <- FOI_eq[i,j]

update(FOI_XL[,,,1]) <- FOI_XL[i,j,k,1] + dt*((lag_rates/dE)*FOI_lag[i,j,k] - (lag_rates/dE)*FOI_XL[i,j,k,1])
update(FOI_XL[,,,2:lag_rates]) <- FOI_XL[i,j,k,l] + dt*((lag_rates/dE)*FOI_XL[i,j,k,l-1] - (lag_rates/dE)*FOI_XL[i,j,k,l])

dim(FOI) <- c(na,nh,num_int)
FOI[,,] <- FOI_XL[i,j,k,lag_rates]


#------------ EIR -----------------
# rate for age group * rate for biting category * FOI for age group * prop of
# infectious mosquitoes
dim(foi_age) <- na
foi_age <- parameter()
dim(rel_foi) <- nh
rel_foi <- parameter()
dim(EIR) <- c(na,nh,num_int)
EIR[,,] <- av_mosq[k] * rel_foi[j] * foi_age[i] * Ivtot/omega

##------------------------------------------------------------------------------
#####################
## IMMUNITY STATES ##
#####################
##------------------------------------------------------------------------------

# See supplementary materials S1 from http://journals.plos.org/plosmedicine/article?id=10.1371/journal.pmed.1000324#s6

# ICM - Maternally acquired immunity acquired by babies from mothers (assumed to be proportional to the immunity of a 15 to 30 year old woman)
# ICA - Immunity acquired due to exposure to previous infection, increases with age
# IC - Clinical immunity. Upon infection, immunity against clinical case. IC = ICA + ICM
# IB - Infection blocking immunity, chances of preventing infection upon receiving infectious bite
# ID - Detection immunity, when immunity suppresses parasite densities this makes it less likely that diagnostics will detect parasite infection

#----------------- DEFINE IMMUNITY PARAMETERS -------------------
dCM <- parameter()
uCA <- parameter()
dCA <- parameter()
dB <- parameter()
uB <- parameter()
dID <- parameter()
uD <- parameter()
x_I <- parameter()
dim(x_I) <- na
age20l <- parameter(type = "integer")
age20u <- parameter(type = "integer")
age_20_factor <- parameter()
PM <- parameter()

#---------- UPDATE IMMUNITY STATES FOR EACH HUMAN COMPARTMENTS ---------

# ICM - maternally acquired immunity
init_ICM <- parameter()
dim(init_ICM) <- c(na,nh,num_int)
initial(ICM[,,]) <- init_ICM[i,j,k]
dim(ICM) <- c(na,nh,num_int)
dim(init_ICM_pre) <- c(nh,num_int)
init_ICM_pre[1:nh,1:num_int] <- PM*(ICA[age20l,i,j] + age_20_factor*(ICA[age20u,i,j]-ICA[age20l,i,j]))

update(ICM[1, 1:nh, 1:num_int]) <- ICM[i,j,k] + dt*(-1/dCM*ICM[i,j,k] + (init_ICM_pre[j,k]-ICM[i,j,k])/x_I[i])
update(ICM[2:na, 1:nh, 1:num_int]) <- ICM[i,j,k] + dt*(-1/dCM*ICM[i,j,k] - (ICM[i,j,k]-ICM[i-1,j,k])/x_I[i])

# ICA - exposure driven immunity
init_ICA <- parameter()
dim(init_ICA) <- c(na,nh,num_int)
initial(ICA[,,]) <- init_ICA[i,j,k]
dim(ICA) <- c(na,nh,num_int)

update(ICA[1, 1:nh, 1:num_int]) <- ICA[i,j,k] + dt*(FOI_smc[i,j,k]/(FOI_smc[i,j,k] * uCA + 1) - 1/dCA*ICA[i,j,k] -ICA[i,j,k]/x_I[i])
update(ICA[2:na, 1:nh, 1:num_int]) <- ICA[i,j,k] + dt*(FOI_smc[i,j,k]/(FOI_smc[i,j,k] * uCA + 1) - 1/dCA*ICA[i,j,k] - (ICA[i,j,k]-ICA[i-1,j,k])/x_I[i])

dim(IC) <- c(na,nh,num_int)
IC[,,] <- ICM[i,j,k] + ICA[i,j,k]

# IB - infection blocking immunity
init_IB <- parameter()
dim(init_IB) <- c(na,nh,num_int)
initial(IB[,,]) <- init_IB[i,j,k]
dim(IB) <- c(na,nh,num_int)

update(IB[1, 1:nh, 1:num_int]) <- IB[i,j,k] + dt*(EIR[i,j,k]/(EIR[i,j,k]* uB + 1) - IB[i,j,k]/dB - IB[i,j,k]/x_I[i])
update(IB[2:na, 1:nh, 1:num_int]) <- IB[i,j,k] + dt*(EIR[i,j,k]/(EIR[i,j,k]* uB + 1) - IB[i,j,k]/dB - (IB[i,j,k]-IB[i-1,j,k])/x_I[i])

# ID - detection immunity
init_ID <- parameter()
dim(init_ID) <- c(na,nh,num_int)
initial(ID[,,]) <- init_ID[i,j,k]
dim(ID) <- c(na,nh,num_int)

update(ID[1, 1:nh, 1:num_int]) <- ID[i,j,k] + dt*(FOI_smc[i,j,k]/(FOI_smc[i,j,k]*uD + 1) - ID[i,j,k]/dID - ID[i,j,k]/x_I[i])
update(ID[2:na, 1:nh, 1:num_int]) <- ID[i,j,k] + dt*(FOI_smc[i,j,k]/(FOI_smc[i,j,k]*uD + 1) - ID[i,j,k]/dID - (ID[i,j,k]-ID[i-1,j,k])/x_I[i])

#----------- PROBABILITIES INFERRED FROM IMMUNITY --------------------
# phi - probability of clinical disease, dependent on clinical immunity
phi0 <- parameter()
phi1 <- parameter()
IC0 <- parameter()
kC <- parameter()
dim(phi) <- c(na,nh,num_int)
phi[1:na,1:nh,1:num_int] <- phi0*((1-phi1)/(1+(IC[i,j,k]/IC0)^kC) + phi1)

# b - probability of disease from infectious bite, depends on infection blocking immunity
b0 <- parameter()
b1 <- parameter()
kB <- parameter()
IB0 <- parameter()
dim(b) <- c(na,nh,num_int)
b[1:na, 1:nh, 1:num_int] <- b0 * ((1-b1)/(1+(IB[i,j,k]/IB0)^kB)+b1)


# p_det - probability of detection by microscopy, immunity decreases chances of
# infection because it pushes parasite density down
aD <- parameter()
fD0 <- parameter()
gammaD <- parameter()
d1 <- parameter()
ID0 <- parameter()
kD <- parameter()
dim(age_vector) <- na
age_vector <- parameter()

dim(fd) <- na
fd[1:na] <- 1-(1-fD0)/(1+(age_vector[i]/aD)^gammaD)
dim(p_det) <- c(na,nh,num_int)
p_det[,,] <- d1 + (1-d1)/(1 + fd[i]*(ID[i,j,k]/ID0)^kD)

##------------------------------------------------------------------------------
#####################
## MOSQUITO STATES ##
#####################
##------------------------------------------------------------------------------

# See supplementary materials S1 from http://journals.plos.org/plosmedicine/article?id=10.1371/journal.pmed.1000324#s6

# Sv - Susceptible mosquitoes
# Pv - Latently infected mosquitoes (Pv is from original version with a single compartment).
# Ev - Latently infected mosquitoes (Ev indicates the explict use of multiple compartments).
# Iv - Infectious mosquitoes

# # Number of pre-infectious compartments (Erlang approximation to sporogony)
# spor_len <- parameter()
#
# # Number of ATN-exposed compartments
# deltaq <- parameter()
# deltaqp1 <- deltaq + 1 # ATN-exposed compartments plus non-exposed

# Number of ATN-exposed compartments (compile-time constants)
deltaq     <- parameter(type = "integer", constant = TRUE)
deltaqp1 <- parameter(type = "integer", constant = TRUE)
spor_len   <- parameter(type = "integer", constant = TRUE)

# Core parameters
gamma_atn  <- parameter()
#xi       <- parameter()  # superseded by Hill kernel for pre-infection blocking
zeta     <- parameter()  # superseded by Hill kernel for EIP suppression
Lambda     <- FOIvdel                # baseline FOI
Lambda00sf <- parameter()            # scale factor for Λ₀⁰
Lambda00   <- Lambda * Lambda00sf    # Λ₀⁰ = Λ * scale factor
rho        <- spor_len / delayMos    # baseline EIP rate
rho_frac   <- parameter()            #
dn0_atn    <- parameter()

# --- Hill kernel parameters (replacing exponential xi / zeta) ---
# Pre-infection blocking (Lambda_i decay across mosquito ATN-exposure compartments)
s_half_pre <- parameter()
nH_pre     <- parameter()

# EIP suppression (rho_i decay across mosquito ATN-exposure compartments)
s_half_eip <- parameter()
nH_eip     <- parameter()

# Post-infection blocking (B_post: probability that an Ev mosquito re-exposed
# to an ATN clears its established infection and returns to Sv[2])
B_max_post  <- parameter()   # asymptotic blocking probability at t_post = 0 (LAB scale)
s_half_post <- parameter()
nH_post     <- parameter()

# --- Bompard TRA → field TBA parameters ---
# m, k from Bompard 2020 / Challenger 2023 SI (defaults for Burkina Faso wild
# mosquitoes: m = 0.000157, k = 4.95e-6).
# Closed form: TBA(tra; m, k) = (b - a) / (1 - a)
#   where a = (k/(k+m))^k, b = (k/(k+m*(1-tra)))^k
m_bompard <- parameter()
k_bompard <- parameter()

# a_bompard depends only on (m, k) so it's constant in time and across compartments
a_bompard <- exp(k_bompard * (log(k_bompard) - log(k_bompard + m_bompard)))

# --- Switches (compile-time integer flags) ---
# use_bompard:  1 = apply Bompard TRA→field TBA, 0 = use lab b directly
# use_eip_hill: 1 = Hill kernel for rho_i, 0 = exponential kernel (legacy form)
use_bompard  <- parameter(type = "logical")#type = "integer", constant = TRUE)
use_eip_hill <- parameter(type = "logical")#type = "integer", constant = TRUE)

# Exp-kernel rate for the EIP suppression (only used when use_eip_hill == 0)
# zeta <- parameter()

# Λ₀(τ) and ρ₀(τ): rho00 is kept here as it depends only on rho.
# In v3, Lambda0_t, rho0_t and dn_atn are no longer single-event analytic
# functions of one t0_atn; they are computed as coverage-weighted averages
# over all ATN distribution events in the "ATN distribution events (v3)"
# block of the INTERVENTION MODEL section below.
rho00 <- rho_frac * rho


# --- Decay of ATN effect across mosquito compartments (eqs. 9–10) ---

# Hill kernel σ(s) = s_half^n / (s_half^n + s^n)
# Lab blocking probability: b_lab(s) = (1 - Λ₀(τ)/Λ) * σ_pre(s)
# Field blocking probability: b_field(s) = Bompard(b_lab(s); m, k)
# Field FOI: Λᵢ(s;τ) = Λ * (1 - b_field(s))
# EIP suppression (no Bompard — sporogony rate, not transmission probability):
#   ρᵢ(s;τ) = ρ - (ρ - ρ₀(τ)) * σ_eip(s)
# where s = i - 0.5  (midpoint of day for compartment i)

dim(Lambda_i)    <- deltaqp1
dim(rho_i)       <- deltaqp1
dim(b_lab_pre)   <- deltaqp1
dim(b_field_pre) <- deltaqp1

# Pre-infection blocking — Hill on lab scale, then (optionally) Bompard → field, then FOI
b_lab_pre[]   <- (1 - Lambda0_t/Lambda) * (s_half_pre^nH_pre / (s_half_pre^nH_pre + (i - 0.5)^nH_pre))
b_field_pre[] <- if (use_bompard) (exp(k_bompard * (log(k_bompard) - log(k_bompard + m_bompard * (1 - b_lab_pre[i])))) - a_bompard) / (1 - a_bompard) else b_lab_pre[i]
Lambda_i[]    <- Lambda * (1 - b_field_pre[i])

# EIP suppression — Hill (use_eip_hill == 1) or exponential (use_eip_hill == 0)
rho_i[] <- if (use_eip_hill) rho - (rho - rho0_t) * (s_half_eip^nH_eip / (s_half_eip^nH_eip + (i - 0.5)^nH_eip)) else rho - (rho - rho0_t) * exp(-zeta * (i - 0.5))

# --- Post-infection blocking (Ev → Sv[2] when re-exposed to ATN) ---
# t_post[j]:      midpoint time-since-infection (in days) for sporogony stage j
# b_lab_post[j]:  lab Hill blocking probability at t_post[j]
# B_post[j]:      field blocking probability after Bompard transformation
# B_post_Ecol[j]: per-stage flux pre-multiplied by B_post[j] for the Sv[2] inflow sum

dim(t_post)       <- spor_len
dim(b_lab_post)   <- spor_len
dim(B_post)       <- spor_len
dim(B_post_Ecol)  <- spor_len

t_post[]      <- (i - 0.5) * delayMos / spor_len
b_lab_post[]  <- B_max_post * (s_half_post^nH_post / (s_half_post^nH_post + t_post[i]^nH_post))
B_post[]      <- if (use_bompard) (exp(k_bompard * (log(k_bompard) - log(k_bompard + m_bompard * (1 - b_lab_post[i])))) - a_bompard) / (1 - a_bompard) else b_lab_post[i]
B_post_Ecol[] <- B_post[i] * Ecol[i]


# Mosquito states
dim(Sv) <- deltaqp1
dim(Ev) <- c(deltaqp1, spor_len)
dim(Iv) <- deltaqp1

#initial state values:
init_Sv <- parameter()
init_Pv <- parameter()
init_Iv <- parameter()
initial(Sv[1]) <- init_Sv * mv0
initial(Sv[2:deltaqp1]) <- 0
#initial(Pv) <- init_Pv * mv0
#initial(Ev[1:spor_len]) <- init_Pv * mv0 / spor_len
#rho <- spor_len / delayMos
kappa <- 1 / deltaq
Ev_ratio <- rho / (rho + mu)
Ev_norm_factor <- (1 - Ev_ratio^spor_len) / (1 - Ev_ratio)
E1_eq <- init_Pv * mv0 / Ev_norm_factor
initial(Ev[1,]) <- E1_eq * Ev_ratio^(j - 1)
initial(Ev[2:deltaqp1,]) <- 0
initial(Iv[1]) <- init_Iv * mv0
initial(Iv[2:deltaqp1]) <- 0

# cA is the infectiousness to mosquitoes of humans in the asmyptomatic compartment broken down
# by age/het/int category, infectiousness depends on p_det which depends on detection immunity
cU <- parameter()
cD <- parameter()
cT <- parameter()
gamma1 <- parameter()
dim(cA) <- c(na,nh,num_int)
cA[,,] <- cU + (cD-cU)*p_det[i,j,k]^gamma1

# Force of infection from humans to mosquitoes
lag_ratesMos <- parameter(type = "integer")

FOIv_eq <- parameter()
dim(FOIv) <- lag_ratesMos
initial(FOIv[]) <- FOIv_eq*delayGam/lag_ratesMos

update(FOIv[1]) <- FOIv[1] + dt*(lag_FOIv - (lag_ratesMos/delayGam)*FOIv[1])
update(FOIv[2:lag_ratesMos]) <- FOIv[i] + dt*((lag_ratesMos/delayGam)*FOIv[i-1] - (lag_ratesMos/delayGam)*FOIv[i])


dim(FOIvijk) <- c(na,nh,num_int)
omega <- parameter()
FOIvijk[1:na, 1:nh, 1:num_int] <- ((cT*smc_rel_c_mask[i,j,k]*T[i,j,k] + cD*smc_rel_c_mask[i,j,k]*D[i,j,k] + cA[i,j,k]*smc_rel_c_mask[i,j,k]*A[i,j,k] + cU*smc_rel_c_mask[i,j,k]*U[i,j,k])/H) *
  rel_foi[j] * av_mosq[k]*foi_age[i]/omega ## For discrete human compartments
lag_FOIv=sum(FOIvijk)

# FOIvdel is equivalent to ince/Sv under previous version
# ince <- FOIv[lag_ratesMos] * lag_ratesMos/delayGam * Sv
FOIvdel <- FOIv[lag_ratesMos] * lag_ratesMos/delayGam

# Updated to account for explicit Erlang EIP

# Current hum->mos FOI depends on the number of individuals now producing gametocytes (12 day lag)
delayGam <- parameter()
delayMos <- parameter()

# Number of mosquitoes born (depends on PL, number of larvae), or is constant outside of seasonality
betaa <- 0.5*PL/dPL

# dim(dSv) <- deltaqp1
# dSv[1] <- dt * (betaa + kappa*Sv[deltaqp1] - (av*delta_atn + (1-delta_atn)*FOIvdel + mu)*Sv[1])
# dSv[2] <- dt * (delta_atn*(av-FOIvdel)*Svtot - (av*delta_atn + (1-delta_atn)*FOIvdel + kappa + mu)*Sv[2])
# dSv[3:deltaqp1] <- dt * (kappa*Sv[i-1] - (av*delta_atn + (1-delta_atn)*FOIvdel + kappa + mu)*Sv[i])

# Susceptible: dS0/dt, dS1/dt, dSi/dt (latex eqs. for S)
dim(dSv) <- deltaqp1
dSv[1] <- dt * ( betaa + kappa*Sv[deltaqp1] - (av*delta_atn + (1 - delta_atn)*Lambda + mu) * Sv[1] )
dSv[2] <- dt * ( delta_atn*(av*(1-dn_atn) - Lambda0_t)*Svtot + av*delta_atn*(1-dn_atn)*sum(B_post_Ecol[1:spor_len]) - (av*delta_atn + (1 - delta_atn)*Lambda_i[2] + kappa + mu) * Sv[2] )
dSv[3:deltaqp1] <- dt * ( kappa*Sv[i-1] - (av*delta_atn + (1 - delta_atn)*Lambda_i[i] + kappa + mu) * Sv[i] )
update(Sv[]) <- if (Sv[i] + dSv[i] < 0) 0 else Sv[i] + dSv[i]


# dim(dEv) <- c(deltaqp1, spor_len)
# dEv[1,1] <- dt * (kappa*Ev[deltaqp1] + (1-delta_atn)*FOIvdel*Sv[1] - (av*delta_atn + rho + mu)*Ev[1,1])
# dEv[2,1] <- dt * (av*delta_atn*Ev + (1-delta_atn)*FOIvdel*Sv[i] - (av*delta_atn + rho + mu)*Ev[i,1])
#
# dEv[,2:spor_len] <- dt * (kappa*Ev[i,j-1] - kappa*Ev[i,j] - mu*Ev[i,j])

dim(Ecol) <- spor_len
Ecol[1:spor_len] <- sum(Ev[,i])
dim(dEv) <- c(deltaqp1, spor_len)

# E_0^1 (i=1, j=1)
dEv[1,1] <- dt * ( kappa*Ev[deltaqp1,1] + (1 - delta_atn)*Lambda*Sv[1] - (av*delta_atn +      rho      + mu) * Ev[1,1] )

# E_1^1 (i=2, j=1)
# Re-exposure inflow scaled by (1 - B_post[1]); the diverted fraction goes to Sv[2].
dEv[2,1] <- dt * ( (1 - B_post[1])*av*delta_atn*(1-dn_atn)*Ecol[1] + delta_atn*(1-dn_atn)*Lambda0_t*Svtot + (1 - delta_atn)*Lambda_i[2]*Sv[2]
                   - (av*delta_atn + kappa + rho_i[2] + mu) * Ev[2,1] )

# E_i^1 (i >= 3, j=1)
dEv[3:deltaqp1,1] <- dt * ( kappa*Ev[i-1,1] + (1 - delta_atn)*Lambda_i[i]*Sv[i]
                            - (av*delta_atn + kappa + rho_i[i] + mu) * Ev[i,1] )

# E_0^j (i=1, j=2..)
dEv[1,2:spor_len] <- dt * ( kappa*Ev[deltaqp1,j] +      rho      *Ev[1,j-1]
                            - (av*delta_atn +      rho      + mu) * Ev[1,j] )

# E_1^j (i=2, j=2..)
# Re-exposure inflow scaled by (1 - B_post[j]); the diverted fraction goes to Sv[2].
# Also adds the (1-dn_atn) survival factor that was missing from the original code,
# matching dEv[2,1] and dIv[2].
dEv[2,2:spor_len] <- dt * ( (1 - B_post[j])*av*delta_atn*(1-dn_atn)*Ecol[j] + rho_i[2]*Ev[2,j-1]
                            - (av*delta_atn + kappa + rho_i[2] + mu) * Ev[2,j] )

# E_i^j (i >= 3, j=2..)
dEv[3:deltaqp1,2:spor_len] <- dt * ( kappa*Ev[i-1,j] + rho_i[i]*Ev[i,j-1]
                                     - (av*delta_atn + kappa + rho_i[i] + mu) * Ev[i,j] )

update(Ev[,]) <- if (Ev[i,j] + dEv[i,j] < 0) 0 else Ev[i,j] + dEv[i,j]


#
# dim(dIv) <- deltaqp1
# dIv[] <- dt * (kappa*Ev[i,spor_len] - mu*Iv[i])
# current total infectious I = sum_i I_i (needed for the i=2 row)
Itot <- sum(Iv[])

dim(dIv) <- deltaqp1
# i = 1  (unexposed row; baseline rho)
dIv[1] <- dt * ( kappa*Iv[deltaqp1] + rho*Ev[1, spor_len] - (av*delta_atn + mu) * Iv[1] )

# i = 2  (exposed 1 day ago; concurrent term αδ I + rho_1 E_1^Δr)
dIv[2] <- dt * ( av*delta_atn*(1-dn_atn)*Itot + rho_i[2]*Ev[2, spor_len]
                 - (av*delta_atn + kappa + mu) * Iv[2] )

# i >= 3  (exposed ≥2 days ago; conveyor κ from previous ATN row + rho_i E_i^Δr)
dIv[3:deltaqp1] <- dt * ( kappa*Iv[i-1] + rho_i[i]*Ev[i, spor_len]
                          - (av*delta_atn + kappa + mu) * Iv[i] )

update(Iv[]) <- if (Iv[i] + dIv[i] < 0) 0 else Iv[i] + dIv[i]

#
# update(Sv[]) <- if(Sv[i] + dSv[i] < 0) 0 else Sv[i] + dSv[i]
# update(Ev[,]) <- if (Ev[i,j] + dEv[i,j] < 0) 0 else Ev[i,j] + dEv[i,j]
# update(Iv[]) <- if(Iv[i] + dIv[i] < 0) 0 else Iv[i] + dIv[i]

# Total mosquito population
initial(Svtot) <- 0
update(Svtot) <- sum(Sv[])

# dim(Evi) <- spor_len
# initial(Evi[]) <- 0
# update(Evi[]) <- sum(Ev[,i])

initial(Evtot) <- 0
update(Evtot) <- sum(Ev[,])

initial(Ivtot) <- 0
update(Ivtot) <- sum(Iv[])

initial(mv) <- 0
update(mv) <- Svtot + Evtot + Ivtot


##------------------------------------------------------------------------------
###################
## LARVAL STATES ##
###################
##------------------------------------------------------------------------------

# Model by White et al.
# (https://parasitesandvectors.biomedcentral.com/articles/10.1186/1756-3305-4-153)

# EL - early larval instar stage
# LL - late larval instar stage
# PL - pupal stage

# mean carrying capacity from initial mosquito density:
dLL <- parameter()
dPL <- parameter()
dEL <- parameter()
muLL <- parameter()
muEL <- parameter()
muPL <- parameter()
gammaL <- parameter()

# fitted entomological parameters:
mv0 <- parameter()
mum <- parameter()
foraging_time <- parameter()
gonotrophic_cycle <- parameter()
betaL <- parameter()

# Entomological variables:
mum_use <- mum
p10 <- exp(-mum_use * foraging_time)  # probability of surviving one feeding cycle
p2 <- exp(-mum_use * gonotrophic_cycle)  # probability of surviving one resting cycle
eov <- betaL/mu*(exp(mu/fv)-1)
beta_larval <- eov*mu*exp(-mu/fv)/(1-exp(-mu/fv)) # Number of eggs laid per day
b_lambda <- (gammaL*muLL/muEL-dEL/dLL+(gammaL-1)*muLL*dEL)
lambda <- -0.5*b_lambda + sqrt(0.25*b_lambda^2 + gammaL*beta_larval*muLL*dEL/(2*muEL*mum_use*dLL*(1+dPL*muPL)))
K0 <- 2*mv0*dLL*mum_use*(1+dPL*muPL)*gammaL*(lambda+1)/(lambda/(muLL*dEL)-1/(muLL*dLL)-1)

# Seasonal carrying capacity KL = base carrying capacity K0 * effect for time of year theta:
KL <- K0 * rain_input
fv <- 1/( foraging_time/(1-zbar) + gonotrophic_cycle ) # mosquito feeding rate (zbar from intervention param.)
mu <- -fv*log(p1*p2) # mosquito death rate

# finding equilibrium and initial values for EL, LL & PL
init_PL <- parameter()
initial(PL) <- init_PL
init_LL <- parameter()
initial(LL) <- init_LL
init_EL <- parameter()
initial(EL) <- init_EL

# (beta_larval (egg rate) * total mosquito) - den. dep. egg mortality - egg hatching
update(EL) <- if(EL + dt*(beta_larval*mv-muEL*(1+(EL+LL)/KL)*EL - EL/dEL) < 0) 0 else (EL + dt*(beta_larval*mv-muEL*(1+(EL+LL)/KL)*EL - EL/dEL))
# egg hatching - den. dep. mortality - maturing larvae
update(LL) <- if(LL + dt*(EL/dEL - muLL*(1+gammaL*(EL + LL)/KL)*LL - LL/dLL) < 0) 0 else (LL + dt*(EL/dEL - muLL*(1+gammaL*(EL + LL)/KL)*LL - LL/dLL))
# pupae - mortality - fully developed pupae
update(PL) <- if(PL + dt*(LL/dLL - muPL*PL - PL/dPL) < 0) 0 else (PL + dt*(LL/dLL - muPL*PL - PL/dPL))


##------------------------------------------------------------------------------
########################
## INTERVENTION MODEL ##
########################
##------------------------------------------------------------------------------

################## SMC #######################
#See supplementary materials of [Thompson, 2022] - https://doi.org/10.1016/S2214-109X(22)00416-8
max_smc_cov <- parameter()

dim(alpha_smc_times, alpha_smc_set) <- parameter(rank = 1)
alpha_smc_times <- parameter(constant = TRUE)
alpha_smc_set <- parameter(constant = TRUE)

alpha_smc <- interpolate(alpha_smc_times, alpha_smc_set, mode = "constant")

dim(alpha_smc_array) <-  c(na,nh,num_int)
alpha_smc_array[,,] <- smc_mask[i,j,k] * alpha_smc #Multiply by proportion of individuals in SMC compartment currently receiving SMC

##Parameters relating to prophylaxis effect of SMC
dim(P_smc_daily) <- n_days + 1
P_smc_daily <- parameter()
P_smc <- interpolate(days, P_smc_daily, "linear")

dim(smc_mask) <- c(na,nh,num_int)
smc_mask <- parameter()
dim(FOI_smc) <- c(na,nh,num_int)
FOI_smc[,,] <- FOI[i,j,k] * (1 - (smc_mask[i,j,k] * P_smc))

##Parameters relating to reduced infectivity of existing infections upon SMC dose
dim(rel_c_days) <- n_days + 1
rel_c_days <- parameter()

rel_c <- interpolate(days, rel_c_days, "constant")
dim(smc_rel_c_mask) <- c(na,nh,num_int)
smc_rel_c_mask[,,] <- 1 - (smc_mask[i,j,k] * (1-rel_c)) #Value = SMC_rel_c for those scheduled to receive SMC on rel_c_days. Else = 1.


################## ITN #######################
# See supplementary materials S2 from http://journals.plos.org/plosmedicine/article?id=10.1371/journal.pmed.1000324#s6
max_itn_cov <- parameter()

Q0 <- parameter()
phi_bednets <- parameter()

# General intervention model terminology:
# r - probability of trying to repeat feed after hitting ITN
# d - probability of dying after hitting ITN
# s - probability of successful feed after hitting ITN

dim(r_itn_daily) <- n_days + 1 #Bednet insecticide decay
dim(s_itn_daily) <- n_days + 1 #Proportion of individuals in the ITN compartment who currently have ITNs.= daily_ITN / max(daily_ITN)
r_itn_daily <- parameter()
s_itn_daily <- parameter()

r_itn <- interpolate(days, r_itn_daily, "linear")
s_itn <- interpolate(days, s_itn_daily, "linear")

################## ATN #######################

# Parameters
# t0_atn     <- parameter()
# Q0_atn     <- parameter()
# lambda_atn <- parameter()
# p_atn      <- parameter()
# phi_atn    <- phi_bednets
#
# # Internal time tracker (odin2 doesn't auto-update 'time')
# initial(tt) <- 0
# update(tt)  <- tt + dt
#
# # ATN coverage state
# initial(Q_atn) <- 0
# #update(Q_atn) <- if (tt < t0_atn) 0 else Q0_atn * exp(-lambda_atn * (tt - t0_atn))
# update(Q_atn) <- if (tt < t0_atn) 0 else if (tt == t0_atn) Q0_atn else Q_atn + dt * (-lambda_atn * Q_atn)
#
#
# # Probability mosquito is exposed to antimalarial drugs on attempted bite
# delta_atn <- p_atn * phi_atn * Q_atn
#
# # Dummy variable so dust2 retains these outputs
# initial(atn_debug) <- 0
# update(atn_debug)  <- delta_atn


# ===================== ATN distribution events (v3) =====================
# v3 generalises the single-event ATN model (v2.1) to an arbitrary number
# of ATN distribution events. t0_atn and Q0_atn are now vectors of length
# n_atn; a single event / scalar usage is simply n_atn = 1 (a length-1
# vector), handled by exactly the same code path. There is no separate
# scalar branch, so the scalar-vs-vector inconsistency present in the ITN
# code (set_bednets.R) is deliberately avoided here.
#
# Mixing rules mirror the ITN treatment in set_bednets.R / utils.R:
#  - Coverage: each new distribution reaches a fraction Q0_atn of the WHOLE
#    population drawn at random (random proportional replacement), so it
#    removes a fraction Q0_atn from every existing category. This is the
#    analogue of get_itn_usage_mat()'s proportional replacement and has a
#    closed form, so no extra state variables are needed.
#  - Drug-effect terms (Lambda0_t, rho0_t, dn_atn): computed per event from
#    that event's own age, then averaged across events weighted by each
#    event's share of current coverage -- the analogue of get_daily_decay()'s
#    usage-weighted averaging of ITN d and r.
#
# Note: lambda_atn, gamma_atn, dn0_atn and p_atn are kept scalar (shared
# across events), matching the request to vectorise only t0_atn and Q0_atn.

# Number of ATN distribution events (compile-time constant; n_atn = 1 for a
# single event / scalar usage).
n_atn <- parameter(type = "integer", constant = TRUE)

# Parameters
dim(t0_atn) <- n_atn
dim(Q0_atn) <- n_atn
t0_atn     <- parameter()   # day of each ATN distribution event (chronological)
Q0_atn     <- parameter()   # initial coverage of each ATN distribution event
lambda_atn <- parameter()   # ATN coverage (retention) decay rate, shared
p_atn      <- parameter()
phi_atn    <- phi_bednets

# --- Per-event coverage with random proportional replacement ---
# repl_factor[i] = product over every later event j that has ALREADY
# occurred of (1 - Q0_atn[j]): the fraction of event i's coverage not yet
# displaced by a subsequent distribution. Computed as exp(sum(log(...)))
# so it works for any n_atn (for n_atn = 1 the sum is empty -> factor 1).
dim(repl_log) <- c(n_atn, n_atn)
repl_log[, ] <- if (t0_atn[j] > t0_atn[i]) (if (t0_atn[j] <= time) log(1 - Q0_atn[j]) else 0) else 0
dim(repl_factor) <- n_atn
repl_factor[] <- exp(sum(repl_log[i, ]))

# Per-event ATN coverage: initial Q0 decayed at lambda_atn since that
# event's t0, then reduced by every later distribution that has since
# occurred. Zero before the event.
dim(Q_atn_each) <- n_atn
Q_atn_each[] <- if (time < t0_atn[i]) 0 else Q0_atn[i] * exp(-lambda_atn * (time - t0_atn[i])) * repl_factor[i]

# Total ATN coverage across all events
Q_atn_t <- sum(Q_atn_each[])

# --- Per-event drug-effect decay (gamma_atn) since each event's t0 ---
dim(Lambda0_each) <- n_atn
dim(rho0_each)    <- n_atn
dim(dn_each)      <- n_atn
Lambda0_each[] <- if (time < t0_atn[i]) Lambda else Lambda - (Lambda - Lambda00) * exp(-gamma_atn * (time - t0_atn[i]))
rho0_each[]    <- if (time < t0_atn[i]) rho    else rho    - (rho    - rho00)    * exp(-gamma_atn * (time - t0_atn[i]))
dn_each[]      <- if (time < t0_atn[i]) 0      else dn0_atn * exp(-gamma_atn * (time - t0_atn[i]))

# --- Coverage-weighted average of the effect terms across events ---
# Weights = each event's share of total current coverage (i.e. conditional
# on a mosquito being exposed to an ATN at all). This is the direct
# analogue of get_daily_decay(): usage-weighted averaging of the per-cohort
# effect. Coverage itself enters separately via delta_atn below, so the
# effect terms are a pure conditional average and are NOT re-scaled by
# coverage. When there is no active coverage they fall back to the no-ATN
# baselines (Lambda, rho, 0).
dim(Lambda0_w) <- n_atn
dim(rho0_w)    <- n_atn
dim(dn_w)      <- n_atn
Lambda0_w[] <- Q_atn_each[i] * Lambda0_each[i]
rho0_w[]    <- Q_atn_each[i] * rho0_each[i]
dn_w[]      <- Q_atn_each[i] * dn_each[i]

Lambda0_t <- if (Q_atn_t > 0) sum(Lambda0_w[]) / Q_atn_t else Lambda
rho0_t    <- if (Q_atn_t > 0) sum(rho0_w[])    / Q_atn_t else rho
dn_atn    <- if (Q_atn_t > 0) sum(dn_w[])      / Q_atn_t else 0

# Probability mosquito is exposed to antimalarial drugs on attempted bite
delta_atn <- p_atn * phi_atn * Q_atn_t

# Track total ATN coverage over time (output like a state variable)
initial(Q_atn) <- 0
update(Q_atn)  <- Q_atn_t

# Dummy variable to ensure delta_atn is kept
initial(atn_debug) <- 0
update(atn_debug)  <- delta_atn
# ========================================================================


################## GENERAL INTERVENTION PARAMETERS #######################
num_int <- parameter()

# cov is a vector of coverages for each intervention category:
dim(cov_) <- 4
cov_[1] <- (1-max_itn_cov)*(1-max_smc_cov)  # {No intervention}
cov_[2] <- max_itn_cov*(1-max_smc_cov) # 	   {ITN only}
cov_[3] <- (1-max_itn_cov)*max_smc_cov	#      {SMC only}
cov_[4] <- max_itn_cov*max_smc_cov #	   {Both itn and SMC}
cov[] <- cov_[i]
dim(cov) <- num_int

# probability that mosquito bites and survives for each intervention category
dim(w_) <- 4
w_[1] <- 1
w_[2] <- 1 - phi_bednets + phi_bednets*s_itn
w_[3] <- 1
w_[4] <- 1 - phi_bednets + phi_bednets*s_itn
w[] <- w_[i]
dim(w) <- num_int

# probability that mosquito is repelled during a single attempt for each int. cat.
dim(z_) <- 4
z_[1] <- 0
z_[2] <- phi_bednets*r_itn
z_[3] <- 0
z_[4] <- phi_bednets*r_itn
z[] <- z_[i]
dim(z) <- num_int

# Calculating Z (zbar) and W (wbar) - see Supplementary materials 2 for details
dim(zhi) <- num_int
dim(whi) <- num_int
zhi[1:num_int] <- cov[i]*z[i]
whi[1:num_int] <- cov[i]*w[i]
zh <- sum(zhi)
wh <- sum(whi)
# Z (zbar) - average probability of mosquito trying again during single feeding attempt
zbar <- Q0*zh
# W (wbar) - average probability of mosquito successfully feeding during single attempt
wbar <- 1 - Q0 + Q0*wh

# p1 is the updated p10 given that interventions are now in place:
p1 <- wbar*p10/(1-zbar*p10)
Q <- 1-(1-Q0)/wbar # updated anthropophagy given interventions
av <- fv*Q # biting rate on humans
dim(av_mosq) <- num_int
av_mosq[1:num_int] <- av*w[i]/wh # rate at which mosquitoes bite each int. cat.


################## GENERAL INTERVENTION PARAMETERS #######################
# Number of intervention combinations:
# 8 when ITN + SMC + ATN are all included
# num_int <- parameter()
#
# # ---------------- COVERAGE CATEGORIES (ITN, SMC, ATN) -----------------
# # Categories:
# # 1: None
# # 2: ITN only
# # 3: SMC only
# # 4: ITN + SMC
# # 5: ATN only
# # 6: ITN + ATN
# # 7: SMC + ATN
# # 8: ITN + SMC + ATN
#
# dim(cov_) <- 8
# cov_[1] <- (1 - max_itn_cov) * (1 - max_smc_cov) * (1 - max_atn_cov)
# cov_[2] <-  max_itn_cov      * (1 - max_smc_cov) * (1 - max_atn_cov)
# cov_[3] <- (1 - max_itn_cov) *  max_smc_cov      * (1 - max_atn_cov)
# cov_[4] <-  max_itn_cov      *  max_smc_cov      * (1 - max_atn_cov)
# cov_[5] <- (1 - max_itn_cov) * (1 - max_smc_cov) *  max_atn_cov
# cov_[6] <-  max_itn_cov      * (1 - max_smc_cov) *  max_atn_cov
# cov_[7] <- (1 - max_itn_cov) *  max_smc_cov      *  max_atn_cov
# cov_[8] <-  max_itn_cov      *  max_smc_cov      *  max_atn_cov
# cov[] <- cov_[i]
# dim(cov) <- num_int
#
# # ---------------- MOSQUITO FEEDING/SURVIVAL EFFECTS -----------------
# # w[i] = probability mosquito bites and survives in category i
# # (ATNs do not currently change feeding survival directly)
# dim(w_) <- 8
# w_[1] <- 1
# w_[2] <- 1 - phi_bednets + phi_bednets * s_itn          # ITN only
# w_[3] <- 1                                              # SMC only
# w_[4] <- 1 - phi_bednets + phi_bednets * s_itn          # ITN + SMC
# w_[5] <- 1                                              # ATN only
# w_[6] <- 1 - phi_bednets + phi_bednets * s_itn          # ITN + ATN
# w_[7] <- 1                                              # SMC + ATN
# w_[8] <- 1 - phi_bednets + phi_bednets * s_itn          # ITN + SMC + ATN
# w[] <- w_[i]
# dim(w) <- num_int
#
# # ---------------- REPELLENCY EFFECTS -----------------
# # z[i] = probability mosquito is repelled during feeding attempt
# # (ATNs assumed to have none)
# dim(z_) <- 8
# z_[1] <- 0
# z_[2] <- phi_bednets * r_itn
# z_[3] <- 0
# z_[4] <- phi_bednets * r_itn
# z_[5] <- 0
# z_[6] <- phi_bednets * r_itn
# z_[7] <- 0
# z_[8] <- phi_bednets * r_itn
# z[] <- z_[i]
# dim(z) <- num_int
#
# # ---------------- WEIGHTED AVERAGES -----------------
# dim(zhi) <- num_int
# dim(whi) <- num_int
# zhi[1:num_int] <- cov[i] * z[i]
# whi[1:num_int] <- cov[i] * w[i]
# zh <- sum(zhi)
# wh <- sum(whi)
#
# # Average probability mosquito retries or succeeds on a feeding attempt
# zbar <- Q0 * zh
# wbar <- 1 - Q0 + Q0 * wh
#
# # Adjusted probabilities with interventions
# p1 <- wbar * p10 / (1 - zbar * p10)
# Q <- 1 - (1 - Q0) / wbar   # updated anthropophagy
# av <- fv * Q               # biting rate on humans
#
# # Category-specific human-biting rate
# dim(av_mosq) <- num_int
# av_mosq[1:num_int] <- av * w[i] / wh


################## ATN COVERAGE TRACKING #######################

# # Parameters from R
# p_atn   <- parameter()
# phi_atn <- parameter()

# Daily ATN coverage vector
# dim(atn_cov_daily) <- n_days + 1
# atn_cov_daily <- parameter()
#
# # Coverage at current time
# Q_atn <- atn_cov_daily[as.integer(time) + 1]

# # Probability mosquito is exposed to antimalarial drug
# delta_atn <- p_atn * phi_atn * Q_atn
#
# # Dummy state so odin2 retains variable
# initial(atn_debug) <- 0
# update(atn_debug)  <- delta_atn



##------------------------------------------------------------------------------
###################
## MODEL OUTPUTS ##
###################
##------------------------------------------------------------------------------
dim(clin_inc) <- c(na,nh,num_int)
initial(clin_inc[,,]) <- init_T[i,j,k]
update(clin_inc[,,]) <- ST_trans[i,j,k] + SD_trans[i,j,k] + AT_trans[i,j,k] + AD_trans[i,j,k] + UD_trans[i,j,k]  + UT_trans[i,j,k]

prev_dim <- parameter()
dim(min_age_prev) <- prev_dim
dim(max_age_prev) <- prev_dim
min_age_prev <- parameter(type = "integer")
max_age_prev <- parameter(type = "integer")

dim(n_prev) <- prev_dim

n_prev[1:prev_dim] <- sum(S[min_age_prev[i]:max_age_prev[i],,]) + sum(T[min_age_prev[i]:max_age_prev[i],,]) + sum(D[min_age_prev[i]:max_age_prev[i],,]) +
  sum(A[min_age_prev[i]:max_age_prev[i],,]) + sum(U[min_age_prev[i]:max_age_prev[i],,]) + sum(P[min_age_prev[i]:max_age_prev[i],,])

dim(detect_prev_full) <- c(na,nh,num_int)
initial(detect_prev_full[,,]) <- init_D[i,j,k]
update(detect_prev_full[,,]) <- T[i,j,k] + D[i,j,k]  + A[i,j,k]*p_det[i,j,k]



dim(detect_prev) <- prev_dim
detect_prev[1:prev_dim] <- sum(detect_prev_full[min_age_prev[i]:max_age_prev[i],,])

##Proportion of population in user defined prevalence age groups
dim(n_ud_prev) <- prev_dim
initial(n_ud_prev[]) <- min_age_prev[i] #user defined prevalence
update(n_ud_prev[]) <- n_prev[i]

##Proportion of population in user defined prevalence age groups with detectable malaria
dim(n_ud_detect_prev) <- prev_dim
initial(n_ud_detect_prev[]) <- min_age_prev[i]
update(n_ud_detect_prev[]) <- detect_prev[i]

# #Incidence
inc_dim <- parameter()
dim(min_age_inc) <- inc_dim
dim(max_age_inc) <- inc_dim
min_age_inc <- parameter(type = "integer")
max_age_inc <- parameter(type = "integer")

dim(n_ud_inc) <- inc_dim #User defined incidence
initial(n_ud_inc[]) <- min_age_inc[i]
update(n_ud_inc[]) <- sum(clin_inc[min_age_inc[i]:max_age_inc[i],,]) / dt

Sh <- sum(S[,,])
Th <- sum(T[,,])
Dh <- sum(D[,,])
Ah <- sum(A[,,])
Uh <- sum(U[,,])
Ph <- sum(P[,,])
H <- Sh + Th + Dh + Ah + Uh + Ph

initial(S_count) <- 0
initial(T_count) <- 0
initial(D_count) <- 0
initial(A_count) <- 0
initial(U_count) <- 0
initial(P_count) <- 0

update(S_count) <- Sh
update(T_count) <- Th
update(D_count) <- Dh
update(A_count) <- Ah
update(U_count) <- Uh
update(P_count) <- Ph

dim(all) <- c(na,nh,num_int)
all[,,] <- S[i,j,k] + T[i,j,k] + D[i,j,k] + A[i,j,k] + U[i,j,k] + P[i,j,k]
##-------------------- IMMUNITY OUTPUTS ----------------------------------
dim(icm_pop) <- c(na,nh,num_int)
icm_pop[,,] <- all[i,j,k] * ICM[i,j,k]
initial(icm_mean) <- 0
update(icm_mean) <- sum(icm_pop[,,]) / H

dim(ica_pop) <- c(na,nh,num_int)
ica_pop[,,] <- all[i,j,k] * ICA[i,j,k]
initial(ica_mean) <- 0
update(ica_mean) <- sum(ica_pop[,,]) / H

dim(id_pop) <- c(na,nh,num_int)
id_pop[,,] <- all[i,j,k] * ID[i,j,k]
initial(id_mean) <- 0
update(id_mean) <- sum(id_pop[,,]) / H

dim(ib_pop) <- c(na,nh,num_int)
ib_pop[,,] <- all[i,j,k] * IB[i,j,k]
initial(ib_mean) <- 0
update(ib_mean) <- sum(ib_pop[,,]) / H

dim(ic_pop) <- c(na,nh,num_int)
ic_pop[,,] <- all[i,j,k] * IC[i,j,k]
initial(ic_mean) <- 0
update(ic_mean) <- sum(ic_pop[,,]) / H

##-------------------- OTHER OUTPUTS ----------------------------------
dim(all_deaths) <- c(na,nh,num_int)
all_deaths[,,] <- S_death[i,j,k] + T_death[i,j,k] + D_death[i,j,k] + A_death[i,j,k] + U_death[i,j,k] + P_death[i,j,k]
initial(natural_deaths) <- 0
update(natural_deaths) <- sum(all_deaths[,,])

dim(epsilon_0) <- c(na,nh,num_int)
epsilon_0[,,] <- (all[i,j,k] * EIR[i,j,k]) / foi_age[i]
initial(EIR_mean) <- 0
update(EIR_mean) <- sum(epsilon_0[,,])

initial(mu_mosq) <- 0
update(mu_mosq) <- mu

##------------- FOR INTERACTION WITH MONTY. ALLOWS FITTING TO PREVALENCE DATA --------------
prevalence <- n_ud_detect_prev[1]/ n_ud_prev[1]
tests <- data()
positive <- data()
positive ~ Binomial(tests, prevalence)

##----------------------- FULL DETECT OUTPUT ----------------------------
dim(detect) <- c(na,nh,num_int)
initial(detect[,,]) <- init_D[i,j,k] + init_T[i,j,k]
update(detect[,,]) <- T[i,j,k] + D[i,j,k]  + A[i,j,k]*p_det[i,j,k]

dim(n) <- c(na,nh,num_int)
initial(n[,,]) <- init_A[i,j,k] + init_T[i,j,k] + init_D[i,j,k] + init_P[i,j,k] + init_U[i,j,k] + init_S[i,j,k]
update(n[,,]) <- all[i,j,k]


