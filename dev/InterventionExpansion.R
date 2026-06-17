# InterventionExpansion.R
# Carry forward the last historical year's intervention rows into the future.
# Returns the site_data list with an extended $interventions data frame.
#
# Arguments:
#   site_data     — a subset site object (from subset_site) with $interventions
#   expand_year   — number of future years to append
#   delay         — not implemented; must be 0
#   counterfactual — not implemented; must be FALSE
expand_interventions <- function(site_data, expand_year, delay = 0, counterfactual = FALSE) {
  if (delay != 0 || isTRUE(counterfactual))
    warning("expand_interventions: delay != 0 / counterfactual != FALSE are not implemented; ignoring.")
  idf      <- site_data$interventions
  last_yr  <- max(idf$year)
  last_rows <- idf[idf$year == last_yr, , drop = FALSE]
  future <- do.call(rbind, lapply(seq_len(expand_year), function(dy) {
    r <- last_rows; r$year <- last_yr + dy; r
  }))
  rownames(future) <- NULL
  site_data$interventions <- rbind(idf, future)
  site_data
}
