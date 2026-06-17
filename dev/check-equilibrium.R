# dev/check-equilibrium.R
# Diagnostic: step the ATN-off adult ODE and compare to expected Erlang equilibrium.
# Run interactively after devtools::load_all().

library(malariasimulation)

# ── Parameters ────────────────────────────────────────────────────────────────
parameters <- get_parameters()           # deltaq=1, spor_len=10, p_atn=0
parameters <- set_equilibrium(parameters, init_EIR = 50.)

deltaq   <- parameters$deltaq
spor_len <- parameters$spor_len
deltaqp1 <- deltaq + 1L
f        <- parameters$blood_meal_rates[[1]]
mu       <- parameters$mum[[1]]

# ── Expected equilibrium ──────────────────────────────────────────────────────
eq <- initial_mosquito_counts(
  parameters, 1L, parameters$init_foim, parameters$total_M
)

# Build compartment names matching the state vector layout
aq_names  <- c("E", "L", "P")
sv_names  <- paste0("Sv[", 0:deltaq, "]")
ev_names  <- as.vector(outer(0:deltaq, 0:(spor_len - 1L),
               FUN = function(q, j) paste0("Ev[", q, ",", j, "]")))
iv_names  <- paste0("Iv[", 0:deltaq, "]")
cpt_names <- c(aq_names, sv_names, ev_names, iv_names)
stopifnot(length(cpt_names) == length(eq))

# ── Build and initialise the solver ──────────────────────────────────────────
timesteps <- 365L
models    <- parameterise_mosquito_models(parameters, timesteps)
solvers   <- parameterise_solvers(models, parameters)

state_init <- solvers[[1]]$get_states()

# ── Step for 1, 10, 100, 365 days ────────────────────────────────────────────
checkpoints <- c(1L, 10L, 100L, 365L)
states_at   <- list()

t <- 0L
for (day in seq_len(timesteps)) {
  kernels <- compute_atn_kernels(day, parameters, parameters$init_foim, 1L)
  adult_mosquito_model_update(
    models[[1]]$.model,
    mu,
    parameters$init_foim,
    f * kernels$delta_atn,
    kernels$delta_atn,
    kernels$dn_atn,
    kernels$Lambda0_t,
    kernels$Lambda_i,
    kernels$rho_i,
    kernels$B_post,
    f
  )
  solvers[[1]]$step()
  t <- t + 1L
  if (t %in% checkpoints) {
    states_at[[as.character(t)]] <- solvers[[1]]$get_states()
  }
}

# ── Print labelled table at each checkpoint ───────────────────────────────────
for (chk in as.character(checkpoints)) {
  actual   <- states_at[[chk]]
  abs_diff <- actual - eq
  rel_diff <- ifelse(abs(eq) > 1e-20, abs_diff / eq, NA_real_)

  tbl <- data.frame(
    compartment = cpt_names,
    expected    = signif(eq,      6),
    actual      = signif(actual,  6),
    abs_diff    = signif(abs_diff, 3),
    rel_diff    = signif(rel_diff, 3),
    stringsAsFactors = FALSE
  )

  # Only print rows where abs relative diff > 1e-8 (suppress exact zeros)
  show <- !is.na(tbl$rel_diff) & abs(tbl$rel_diff) > 1e-8
  if (!any(show)) {
    cat(sprintf("\n=== Day %s: all compartments within 1e-8 ===\n", chk))
  } else {
    cat(sprintf(
      "\n=== Day %s: %d compartments with |rel_diff| > 1e-8 ===\n",
      chk, sum(show)
    ))
    print(tbl[show, ], row.names = FALSE)
  }
}

# ── Summary: max relative drift over all compartments at day 365 ─────────────
actual_365 <- states_at[["365"]]
rel_365    <- abs((actual_365 - eq) / ifelse(abs(eq) > 1e-20, eq, 1))
cat(sprintf(
  "\nMax |rel diff| at day 365 : %.3e  (compartment: %s)\n",
  max(rel_365, na.rm = TRUE),
  cpt_names[which.max(rel_365)]
))
cat(sprintf(
  "Max |abs diff| at day 365 : %.3e  (compartment: %s)\n",
  max(abs(actual_365 - eq), na.rm = TRUE),
  cpt_names[which.max(abs(actual_365 - eq))]
))

# ── Check: does the initial state match eq exactly? ───────────────────────────
cat("\n=== Init state vs equilibrium ===\n")
init_rel <- abs((state_init - eq) / ifelse(abs(eq) > 1e-20, eq, 1))
cat(sprintf("Max |rel diff| init vs eq  : %.3e\n", max(init_rel, na.rm = TRUE)))
if (max(init_rel, na.rm = TRUE) > 1e-10) {
  show_init <- init_rel > 1e-10
  cat("Compartments with init mismatch:\n")
  print(data.frame(
    compartment = cpt_names[show_init],
    expected    = eq[show_init],
    actual_init = state_init[show_init],
    rel_diff    = init_rel[show_init]
  ), row.names = FALSE)
}
