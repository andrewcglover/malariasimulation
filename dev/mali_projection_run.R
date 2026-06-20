# =====================================================================
# mali_projection_run.R
#   Mali admin-1 (no urban/rural split): 6-year forward projections under
#   four future net arms (none/cfp/atn/pyr_atn), median single run, with
#   continuous distribution (CD) of future nets. Parallel on Windows.
#
#   Past = calibrated site history (kept). Future non-net interventions are
#   controlled by the `future_interventions` toggle list near line 35; defaults
#   carry CM + SMC forward and zero IRS / vaccines / PMC / LSM. Nets built
#   manually (past from site df + future CD per arm) via one set_bednets().
#
#   Run from the FORK ROOT after building the fork interactively:
#       devtools::load_all()        # compiles the .dll once
#   then source this file.
# =====================================================================

library(malariasimulation)
library(site)
library(dplyr); library(tidyr)
library(parallel)
source("dev/InterventionExpansion.R")

# ---------------------------------------------------------------------
# 0. Paths & top-level settings
# ---------------------------------------------------------------------
SITE_FILE <- "dev/site_files/without_split/MLI.rds"
OUT_FILE  <- "dev/outputs/mali_projection_results.rds"
FORK_PATH <- normalizePath(".")          # for pkgload::load_all() in workers
N_CORES   <- min(16L, max(1L, parallel::detectCores() - 1L))

human_pop      <- 100000L  # per-region population (smooth single-run incidence)
n_future_years <- 6L
arms           <- c("none", "cfp", "atn", "pyr_atn")
retention_override <- NULL # numeric (days) to override site mean_retention (default: use site)
deltaq_use     <- 10L

# Future non-net site interventions: carried forward at 2024 levels if TRUE,
# zeroed over the future window if FALSE. Historical rows are ALWAYS kept.
# Net columns (itn_*) are always zeroed here — nets are rebuilt per-arm below.
# DEFAULT reproduces the current behaviour exactly (CM + SMC on, rest off).
# To run with future CM/SMC off, set case_management and smc to FALSE.
future_interventions <- list(
  case_management = TRUE,   # tx_cov
  smc             = TRUE,   # smc_cov
  irs             = FALSE,  # irs_cov
  rtss            = FALSE,  # rtss_cov
  r21             = FALSE,  # r21_cov
  pmc             = FALSE,  # pmc_cov
  lsm             = FALSE   # lsm_cov
)

# Load the site early: net retention is sourced from it (used by the CD math below).
site_obj <- readRDS(SITE_FILE)

# Net retention (days): site mean_retention unless manually overridden.
# NB: site mean_retention (~2014 d) is much longer than the previous hardcoded
# 588 d and materially changes the CD top-up coverage below — intended, site-driven.
site_retention <- unique(site_obj$interventions$mean_retention)
stopifnot(length(site_retention) == 1L)
retention_time <- if (!is.null(retention_override)) retention_override else site_retention
message(sprintf("Net retention = %.1f days (%s)", retention_time,
                if (is.null(retention_override)) "site mean_retention" else "manual override"))

# ---------------------------------------------------------------------
# 1. CD coverage math (identical to atn_local_test_cd_v1.R)
# ---------------------------------------------------------------------
campaign_cov <- 0.9          # campaign COVERAGE (fraction receiving a new net)
CD_frac      <- 0.155        # CD share of post-campaign net use
cd_interval  <- 365 / 12    # monthly top-up cadence

d_cd     <- exp(-cd_interval / retention_time)
cd_floor <- CD_frac * campaign_cov / (1 - CD_frac * (1 - campaign_cov))
cd_cov   <- cd_floor * (1 - d_cd) / (1 - cd_floor * d_cd)
message(sprintf("CD: campaign_cov=%.2f -> floor=%.4f, monthly_top-up_cov=%.4f",
                campaign_cov, cd_floor, cd_cov))

# ---------------------------------------------------------------------
# 2. Median net-efficacy params
#    gamman in the RDS files is in years; multiply by 365 -> days.
# ---------------------------------------------------------------------
cfp_pars  <- readRDS("dev/itn_params/dat_res_cfp.rds")  |> mutate(gamman = gamman * 365)
only_pars <- readRDS("dev/itn_params/dat_res_only.rds") |> mutate(gamman = gamman * 365)

# Pick net params at the nearest resistance-grid level (grid at 0.01 increments).
med_net <- function(pars, res) {
  grid    <- sort(unique(pars$resistance))
  res_use <- grid[which.min(abs(grid - res))]
  pars |> filter(resistance == res_use) |>
    summarise(dn0 = median(dn0), rn0 = median(rn0), gamman = median(gamman))
}

# ---------------------------------------------------------------------
# 3. Median ATN kernel parameters (Stan posteriors)
# ---------------------------------------------------------------------
eip_samples <- rstan::extract(readRDS("dev/atn_params/eip_fit_hill_result.rds"))
tra_samples <- rstan::extract(readRDS("dev/atn_params/atn_fit_tra_main_bidir_splitnH.rds"))
atn_kern <- list(
  s_half_eip  = median(eip_samples$s_half),
  nH_eip      = median(eip_samples$nH),
  rho_frac    = median(eip_samples$r0_frac),
  s_half_pre  = median(tra_samples$s_half_pre),
  nH_pre      = median(tra_samples$nH_pre),
  B_max_post  = median(tra_samples$b_max),
  s_half_post = median(tra_samples$s_half_post),
  nH_post     = median(tra_samples$nH_post)
)
atn_halflife <- median(c(only_pars$gamman, cfp_pars$gamman))  # days
gamma_atn    <- log(2) / atn_halflife

render_overrides <- list(
  human_population                      = human_pop,
  prevalence_rendering_min_ages         = 2 * 365,
  prevalence_rendering_max_ages         = 10 * 365,
  age_group_rendering_min_ages          = 2 * 365,
  age_group_rendering_max_ages          = 10 * 365,
  clinical_incidence_rendering_min_ages = 0,
  clinical_incidence_rendering_max_ages = 100 * 365
)
form_overrides <- list(use_eip_hill = TRUE, use_bompard = TRUE)

# ---------------------------------------------------------------------
# 4. Load site; establish the calendar
#    baseline_year = min(year), step_t = itn_dist_day + (year - baseline_year)*365
#    (confirmed from site:::add_time and site:::add_itns source)
# ---------------------------------------------------------------------
# site_obj already loaded above (for site-sourced net retention).

# Canonical region list from the sites table (not interventions, to be safe).
regions    <- site_obj$sites$name_1                 # 9 admin-1 regions
start_year <- min(site_obj$interventions$year)       # 2000
hist_last  <- max(site_obj$interventions$year)       # 2024
future_yr0 <- hist_last + 1L                         # 2025
n_years    <- (hist_last - start_year + 1L) + n_future_years   # 31
n_steps    <- n_years * 365L                         # 11315
future_start_day <- (future_yr0 - start_year) * 365L # 9125

# Future mass-campaign days: every 36 months (1095 days) within the 6-yr window.
# 36 months * (365/12 days/month) = 1095 days exactly.
future_campaign_days <- seq(future_start_day,
                            future_start_day + (n_future_years - 1L) * 365L,
                            by = 3L * 365L)

# ---------------------------------------------------------------------
# 5. Net schedule builders
# ---------------------------------------------------------------------

# PAST: historical ITN distributions from the site df (rows with itn_input_dist > 0).
# gamman in the site df is in years; * 365 -> days (same conversion as add_itns).
build_past_schedule <- function(idf) {
  pd <- idf[!is.na(idf$itn_input_dist) & idf$itn_input_dist > 0, ]
  if (nrow(pd) == 0L) return(NULL)
  list(
    timesteps = (pd$year - start_year) * 365L + pd$itn_distribution_day,
    coverages = pd$itn_input_dist,
    dn0 = pd$dn0,
    rn  = pd$rn0,   # site df column is rn0; set_bednets param is rn
    rnm = pd$rnm,
    gam = pd$gamman * 365   # site df gamman is in years
  )
}

# FUTURE: monthly CD grid + 3-yearly campaigns, arm-specific net type.
# Net efficacy (dn0/rn/gamman) is projected per distribution year using the site
# file's year-by-year pyrethroid resistance (site_obj$vectors$pyrethroid_resistance,
# 2000-2050), taking the median posterior draw conditional on that year's resistance.
build_future_schedule <- function(arm, region) {
  if (arm == "none") return(NULL)   # no future nets for the 'none' arm

  grid    <- unique(round(seq(future_start_day, n_steps, by = cd_interval)))
  is_camp <- vapply(grid, function(t) any(abs(t - future_campaign_days) < 1L), logical(1))
  cov     <- ifelse(is_camp, campaign_cov, cd_cov)
  n       <- length(grid)

  # Projected resistance at each distribution's calendar year, for this region.
  rtab      <- site_obj$vectors$pyrethroid_resistance
  rtab      <- rtab[rtab$name_1 == region, ]
  grid_year <- start_year + (grid %/% 365L)
  res_grid  <- rtab$pyrethroid_resistance[match(grid_year, rtab$year)]

  if (arm == "cfp") {                     # pyrethroid + CFP net
    pars <- lapply(res_grid, function(r) med_net(cfp_pars, r))
    dn0  <- vapply(pars, `[[`, numeric(1), "dn0")
    rn   <- vapply(pars, `[[`, numeric(1), "rn0")
    rnm  <- rep(0.24, n)
    gam  <- vapply(pars, `[[`, numeric(1), "gamman")
  } else if (arm == "atn") {              # antimalarial net, no insecticide
    dn0 <- rep(0, n); rn <- rep(0.24, n); rnm <- rep(0.24 - 1e-9, n)
    gam <- rep(atn_halflife, n)
  } else if (arm == "pyr_atn") {          # pyrethroid + antimalarial
    pars <- lapply(res_grid, function(r) med_net(only_pars, r))
    dn0  <- vapply(pars, `[[`, numeric(1), "dn0")
    rn   <- vapply(pars, `[[`, numeric(1), "rn0")
    rnm  <- rep(0.24, n)
    gam  <- vapply(pars, `[[`, numeric(1), "gamman")
  }
  list(timesteps = grid, coverages = cov,
       dn0 = dn0, rn = rn, rnm = rnm, gam = gam)
}

# ---------------------------------------------------------------------
# 6. Build parameters for one region x arm
# ---------------------------------------------------------------------
build_params <- function(region, arm) {

  # Extract the one admin-1 row from the sites table and subset the site object.
  # subset_site() is the correct extractor (site package exports subset_site, not single_site).
  site_row <- site_obj$sites[site_obj$sites$name_1 == region, , drop = FALSE]
  ms       <- site::subset_site(site_obj, site_row)

  # 6a. Extend future: zero net columns (rebuilt per-arm); apply `future_interventions`
  #     toggles for non-net columns.  Historical rows are always kept unchanged.
  ms_ext <- expand_interventions(ms, expand_year = n_future_years)
  fut    <- ms_ext$interventions$year >= future_yr0

  # Net columns — always zeroed (rebuilt from scratch by the arm-specific set_bednets call).
  net_zero_cols <- c("itn_input_dist", "itn_use")
  for (col in net_zero_cols)
    if (col %in% names(ms_ext$interventions))
      ms_ext$interventions[[col]][fut] <- 0

  # Toggleable non-net interventions (zero future window when toggled FALSE).
  toggle_cols <- list(
    case_management = "tx_cov",  smc  = "smc_cov", irs = "irs_cov",
    rtss            = "rtss_cov", r21 = "r21_cov", pmc = "pmc_cov",
    lsm             = "lsm_cov"
  )
  for (nm in names(toggle_cols)) {
    if (!isTRUE(future_interventions[[nm]])) {
      col <- toggle_cols[[nm]]
      if (col %in% names(ms_ext$interventions))
        ms_ext$interventions[[col]][fut] <- 0
    }
  }

  # 6b. ATN construction-time overrides (antimalarial arms only).
  #     These pass through site_parameters -> get_parameters(overrides=...).
  #     (Confirmed from site:::site_parameters source.)
  fsch <- build_future_schedule(arm, region)
  atn_overrides <- list()
  if (arm %in% c("atn", "pyr_atn")) {
    atn_overrides <- c(list(
      p_atn     = 0.9,
      deltaq    = deltaq_use,
      gamma_atn = gamma_atn,
      Q0_atn    = fsch$coverages,
      t0_atn    = fsch$timesteps,
      n_atn     = length(fsch$timesteps)
    ), atn_kern)
    # Pyrethroid mortality for pyr_atn flows through the ITN-side dn0/rn in the
    # net schedule (resistance-projected); the ATN kernel's dn0_atn represents
    # only the antimalarial's extra mortality (default 0).
  }

  # 6c. Build calibrated parameters from the extended site unit.
  #
  #   CONFIRMED API (from source inspection):
  #   - vectors   = ms$vectors$vector_species  (flat species df, not the whole list)
  #   - seasonality = ms$seasonality$seasonality_parameters  (Fourier df, not the list)
  #   - eir       = ms$eir$eir  (numeric scalar)
  #   - overrides are forwarded directly to get_parameters(overrides=...)
  #   - set_bednets is called internally by add_itns for historical rows;
  #     our manual call below OVERWRITES it completely.
  #
  p <- site::site_parameters(
    interventions = ms_ext$interventions,
    demography    = ms_ext$demography,
    vectors       = ms_ext$vectors$vector_species,
    seasonality   = ms_ext$seasonality$seasonality_parameters,
    eir           = ms$eir$eir,
    overrides     = c(render_overrides, form_overrides, atn_overrides,
                      list(ode_max_steps = 1e8,
                           a_tol = 0.1))   # raised from 1e-4: prevents near-zero thrashing
                                           # on Erlang adult-mosquito ODE at low density
                                           # (none/cfp: 27 states, deltaq=1; atn/pyr_atn: 135).
                                           # Without the C++ non-negativity floor, a_tol=0.1
                                           # caused NaN cascades (negative overshoot in dry
                                           # season trough). Floor + a_tol=0.1 together fix it.
                                           # Max Ivtot error ~0.1–1 during near-elimination
                                           # (biologically inert; human IBM reseeds when EIR
                                           # recovers). See dev/bamako_solver_diagnosis.md.
  )

  # 6d. Combined net schedule: past (from site df) + future (CD, arm-specific).
  #     Concatenate, sort chronologically, then pass as a single set_bednets call.
  #     This OVERWRITES whatever set_bednets the site_parameters call set.
  past <- build_past_schedule(ms$interventions)

  ts  <- c(past$timesteps,  fsch$timesteps)
  cov <- c(past$coverages,  fsch$coverages)
  dn0 <- c(past$dn0,        fsch$dn0)
  rn  <- c(past$rn,         fsch$rn)
  rnm <- c(past$rnm,        fsch$rnm)
  gam <- c(past$gam,        fsch$gam)
  # Note: c(x, NULL) = x in R, so the 'none' arm (fsch=NULL) reduces to past only.

  ord <- order(ts)
  ts <- ts[ord]; cov <- cov[ord]; dn0 <- dn0[ord]
  rn <- rn[ord]; rnm <- rnm[ord]; gam <- gam[ord]

  # set_bednets requires dn0/rn/rnm as matrices with ncol = number of species.
  # (Confirmed from set_bednets source: "rows need to align with species".)
  n_sp <- length(p$species)   # 3 for MLI (gambiae / arabiensis / funestus)
  p <- set_bednets(p,
    timesteps = ts,
    coverages = cov,
    retention = retention_time,
    dn0       = matrix(rep(dn0, n_sp), ncol = n_sp),
    rn        = matrix(rep(rn,  n_sp), ncol = n_sp),
    rnm       = matrix(rep(rnm, n_sp), ncol = n_sp),
    gamman    = gam
  )
  p
}

# ---------------------------------------------------------------------
# 7. Run one region x arm row
# ---------------------------------------------------------------------
run_one <- function(row) {
  region <- row$region; arm <- row$arm
  tryCatch({
    p <- build_params(region, arm)
    r <- run_simulation(timesteps = n_steps, parameters = p)
    as.data.frame(r) |>
      mutate(
        region    = region,
        arm       = arm,
        year_rel  = (timestep - future_start_day) / 365,
        pfpr2to10 = n_detect_lm_730_3649 / n_age_730_3649,
        clin_inc  = n_inc_clinical_0_1824 + n_inc_clinical_1825_5474 + n_inc_clinical_5475_36499,
        EIR_gambiae_pp    = EIR_gambiae    / human_pop,
        EIR_arabiensis_pp = EIR_arabiensis / human_pop,
        EIR_funestus_pp   = EIR_funestus   / human_pop
      )
  }, error = function(e) {
    list(.__error__ = TRUE, region = region, arm = arm, message = conditionMessage(e))
  })
}

# ---------------------------------------------------------------------
# 8. Parallel sweep (Windows PSOCK; workers reuse the .dll built interactively)
#    Skipped when options(mali_test_mode = TRUE) — lets test scripts source
#    this file to get functions/data without launching the cluster.
# ---------------------------------------------------------------------
if (isTRUE(getOption("mali_test_mode"))) {
  message("mali_test_mode = TRUE: setup complete, skipping parallel sweep.")
} else {

grid_df <- expand.grid(region = regions, arm = arms, stringsAsFactors = FALSE)
rows    <- split(grid_df, seq_len(nrow(grid_df)))

cl <- parallel::makeCluster(N_CORES)
on.exit(parallel::stopCluster(cl), add = TRUE)

parallel::clusterExport(cl, "FORK_PATH")
parallel::clusterEvalQ(cl, {
  suppressMessages({ library(site); library(dplyr); library(tidyr) })
  pkgload::load_all(FORK_PATH, quiet = TRUE, compile = FALSE)
})
parallel::clusterExport(cl, c(
  "build_params", "build_past_schedule", "build_future_schedule", "med_net",
  "expand_interventions", "run_one",
  "site_obj", "regions", "start_year", "hist_last", "future_yr0",
  "future_start_day", "future_campaign_days", "n_steps", "n_future_years",
  "cfp_pars", "only_pars", "atn_kern", "atn_halflife", "gamma_atn",
  "cd_cov", "campaign_cov", "cd_interval", "retention_time", "deltaq_use",
  "render_overrides", "form_overrides", "human_pop",
  "future_interventions"
))

t0 <- proc.time()[["elapsed"]]
res_list <- parallel::parLapplyLB(cl, rows, run_one)
message(sprintf("Done %d runs in %.0f s", length(rows),
                proc.time()[["elapsed"]] - t0))

is_error <- vapply(res_list, function(x) isTRUE(x$.__error__), logical(1))
if (any(is_error)) {
  failures <- dplyr::bind_rows(lapply(res_list[is_error], function(x)
    data.frame(region = x$region, arm = x$arm, message = x$message,
               stringsAsFactors = FALSE)))
  message(sprintf("WARNING: %d job(s) failed and excluded from output:", sum(is_error)))
  print(failures)
  write.csv(failures,
            sub("\\.rds$", "_failures.csv", OUT_FILE),
            row.names = FALSE)
}
df_full <- dplyr::bind_rows(res_list[!is_error])
dir.create(dirname(OUT_FILE), recursive = TRUE, showWarnings = FALSE)
saveRDS(list(
  results        = df_full,
  shape          = site_obj$shape$level_1,   # sf polygon object for admin-1
  meta           = list(
    start_year     = start_year,
    hist_last      = hist_last,
    future_yr0     = future_yr0,
    future_start_day = future_start_day,
    n_future_years = n_future_years,
    human_pop      = human_pop,
    cd_floor       = cd_floor,
    cd_cov         = cd_cov
  )
), OUT_FILE)
message("Saved -> ", OUT_FILE)

} # end if (!mali_test_mode)
