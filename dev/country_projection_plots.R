# =====================================================================
# country_projection_plots.R
#   Reads the results rds produced by country_projection_run.R and produces:
#     (1) prevalence-over-time, geographically facetted by region (geofacet)
#     (2) three clinical-cases-averted choropleths:
#           Pyr-CFP vs none, ATN vs none, ATN additional vs Pyr-CFP
#     (2b) resistance x baseline-EIR scatter, coloured by ATN - Pyr-CFP
#   Change COUNTRY_ISO to match the country_projection_run.R setting.
#   Run from the fork root after country_projection_run.R has completed.
# =====================================================================

library(dplyr); library(tidyr); library(ggplot2); library(zoo)
library(sf); library(patchwork)

# ---------------------------------------------------------------------
# 0. Config — must match country_projection_run.R
# ---------------------------------------------------------------------
COUNTRY_ISO  <- "BFA"
RESULTS_FILE <- sprintf("dev/outputs/%s_projection_results.rds", tolower(COUNTRY_ISO))
SITE_FILE    <- sprintf("dev/site_files/without_split/%s.rds", COUNTRY_ISO)

obj  <- readRDS(RESULTS_FILE)
df   <- obj$results
shape <- obj$shape
meta <- obj$meta

country_name <- if (!is.null(meta$country_name)) meta$country_name else COUNTRY_ISO

# Shape column holding region names.
SHAPE_KEY <- "name_1"

arm_labels <- c(none    = "No future nets",
                cfp     = "Future Pyr-CFP",
                atn     = "Future ATN",
                pyr_atn = "Future Pyr-ATN")
arm_cols   <- c("No future nets" = "grey50",
                "Future Pyr-CFP" = "#009988",
                "Future ATN"     = "#EE7733",
                "Future Pyr-ATN" = "#CC3311")

# Dynamic excluded-region note (replaces hardcoded "Bamako excluded" in the Mali script).
missing_regions <- setdiff(shape[[SHAPE_KEY]], unique(df$region))
excl_note <- if (length(missing_regions) > 0L)
  paste0(" Excluded (no successful run): ", paste(missing_regions, collapse = ", "), ".")
else
  ""

# ---------------------------------------------------------------------
# 1. Geographic facets of prevalence over the future window
# ---------------------------------------------------------------------
df_win <- df |>
  filter(year_rel >= 0, year_rel <= meta$n_future_years) |>
  mutate(arm_f = factor(arm_labels[arm], levels = arm_labels))

country_grid_raw <- geofacet::grid_auto(shape, names = SHAPE_KEY, seed = 1)
names(country_grid_raw)[names(country_grid_raw) == paste0("name_", SHAPE_KEY)] <- "name"
country_grid_raw$code <- country_grid_raw$name
present_regions <- unique(df_win$region)
country_grid <- country_grid_raw[country_grid_raw$code %in% present_regions,
                                 c("row", "col", "code", "name")]

df_plot <- df |>
  filter(year_rel >= -3) |>
  mutate(
    arm_f    = factor(arm_labels[arm], levels = arm_labels),
    cal_year = year_rel + meta$future_yr0,
    clin_rate = clin_inc / meta$human_pop * 1000
  ) |>
  group_by(region, arm_f) |>
  arrange(cal_year, .by_group = TRUE) |>
  mutate(
    pfpr_roll      = zoo::rollmean(pfpr2to10 * 100, k = 365L, fill = NA, align = "center"),
    clin_rate_roll = zoo::rollmean(clin_rate,        k = 365L, fill = NA, align = "center")
  ) |>
  ungroup()

vline_df    <- data.frame(xintercept = meta$future_yr0)
plot_caption <- paste0("Dashed line = future distribution start (", meta$future_yr0, ").",
                       " Thick line = 365-day rolling mean.",
                       "\n'No future nets': historical ITNs only, decaying from ", meta$future_yr0,
                       " (no new distributions).", excl_note)

p_clin_series <- ggplot(df_plot, aes(cal_year, clin_rate, colour = arm_f)) +
  geom_vline(data = vline_df, aes(xintercept = xintercept),
             linetype = "dashed", colour = "grey40", linewidth = 0.4) +
  geom_line(linewidth = 0.4, alpha = 0.2) +
  geom_line(aes(y = clin_rate_roll), linewidth = 0.8, na.rm = TRUE) +
  scale_colour_manual(values = arm_cols) +
  geofacet::facet_geo(~ region, grid = country_grid) +
  theme_minimal(base_size = 11) +
  labs(x = "Year",
       y = "Clinical incidence (per 1,000 / day, all ages)",
       colour = "",
       title = paste0(country_name, " — projected clinical incidence by region and future net scenario"),
       caption = plot_caption)

p_prev <- ggplot(df_plot, aes(cal_year, pfpr2to10 * 100, colour = arm_f)) +
  geom_vline(data = vline_df, aes(xintercept = xintercept),
             linetype = "dashed", colour = "grey40", linewidth = 0.4) +
  geom_line(linewidth = 0.4, alpha = 0.2) +
  geom_line(aes(y = pfpr_roll), linewidth = 0.8, na.rm = TRUE) +
  scale_colour_manual(values = arm_cols) +
  geofacet::facet_geo(~ region, grid = country_grid) +
  theme_minimal(base_size = 11) +
  labs(x = "Year",
       y = expression(italic(Pf) * PR[2 - 10] * " (%)"),
       colour = "",
       title = paste0(country_name, " — projected prevalence by region and future net scenario"),
       caption = plot_caption)

out_prev  <- sprintf("dev/outputs/%s_prevalence_facet.png",  tolower(COUNTRY_ISO))
out_clin  <- sprintf("dev/outputs/%s_incidence_facet.png",   tolower(COUNTRY_ISO))
ggsave(out_prev, p_prev,         width = 11, height = 8, dpi = 150)
ggsave(out_clin, p_clin_series,  width = 11, height = 8, dpi = 150)

# ---------------------------------------------------------------------
# 2. Clinical cases averted (future window) per 1,000 pop per year
# ---------------------------------------------------------------------
per1000_yr <- function(total) total / meta$human_pop * 1000 / meta$n_future_years

tot <- df_win |>
  group_by(region, arm) |>
  summarise(cases = sum(clin_inc), .groups = "drop") |>
  pivot_wider(names_from = arm, values_from = cases) |>
  mutate(
    cfp_averted    = per1000_yr(none - cfp),
    atn_averted    = per1000_yr(none - atn),
    atn_additional = per1000_yr(cfp  - atn)
  )

map_df <- shape |> left_join(tot, by = setNames("region", SHAPE_KEY))

make_map <- function(fill_col, title, diverging = FALSE) {
  g <- ggplot(map_df) +
    geom_sf(aes(fill = .data[[fill_col]]), colour = "white", linewidth = 0.2) +
    theme_void(base_size = 11) +
    labs(title = title, fill = "Averted /\n1000 / yr")
  if (diverging) {
    g + scale_fill_gradient2(low = "#B2182B", mid = "grey95", high = "#2166AC",
                             midpoint = 0)
  } else {
    g + scale_fill_viridis_c(option = "D", direction = 1)
  }
}

rng   <- range(c(tot$cfp_averted, tot$atn_averted), na.rm = TRUE)
m_cfp <- make_map("cfp_averted",    "Pyr-CFP vs\nno future nets") +
           scale_fill_viridis_c(limits = rng, option = "D")
m_atn <- make_map("atn_averted",    "ATN vs\nno future nets") +
           scale_fill_viridis_c(limits = rng, option = "D")
m_add <- make_map("atn_additional", "ATN additional\n(vs Pyr-CFP)", diverging = TRUE)

maps_caption <- paste0(
  "'No future nets' = historical ITNs decaying from ", meta$future_yr0,
  " with no replacements; incidence rebounds toward unprotected equilibrium.\n",
  "High averted counts reflect this rebound, not the current protected baseline.",
  excl_note
)

p_maps <- (m_cfp | m_atn | m_add) +
  patchwork::plot_annotation(
    title   = paste0(country_name,
                     " — clinical cases averted over ", meta$n_future_years,
                     "-year projection (per 1,000 pop / yr)"),
    caption = maps_caption
  )

out_maps <- sprintf("dev/outputs/%s_averted_maps.png", tolower(COUNTRY_ISO))
ggsave(out_maps, p_maps, width = 14, height = 5, dpi = 150)

message(sprintf("Saved: %s, %s, %s", out_prev, out_clin, out_maps))

# ---------------------------------------------------------------------
# 2b. Resistance x baseline-EIR scatter, coloured by ATN - Pyr-CFP
#     change in annual clinical incidence (per 1,000 / yr).
# ---------------------------------------------------------------------
site_obj2 <- readRDS(SITE_FILE)

fut_years <- meta$future_yr0 + 0:(meta$n_future_years - 1L)
res_fut <- site_obj2$vectors$pyrethroid_resistance |>
  filter(year %in% fut_years) |>
  group_by(name_1) |>
  summarise(res_fut = mean(pyrethroid_resistance), .groups = "drop")

eir_base <- site_obj2$eir |>
  filter(sp == "pf") |>
  select(name_1, eir)

atn_vs_cfp <- df_win |>
  filter(arm %in% c("atn", "cfp")) |>
  group_by(region, arm) |>
  summarise(cases = sum(clin_inc), .groups = "drop") |>
  pivot_wider(names_from = arm, values_from = cases) |>
  mutate(delta_clin = per1000_yr(atn - cfp),
         pct_clin   = (atn - cfp) / cfp * 100)

scatter_df <- atn_vs_cfp |>
  left_join(eir_base, by = c(region = "name_1")) |>
  left_join(res_fut,  by = c(region = "name_1"))

fill_lim <- max(abs(scatter_df$delta_clin), na.rm = TRUE) * c(-1, 1)

# Manual limit for the % change scatter (single positive number, e.g. 50 for ±50%).
# NULL = auto (symmetric at max observed abs value). Out-of-range points are
# squished to the palette extremes.
pct_scatter_lim <- 40 # NULL

if (is.null(pct_scatter_lim)) {
  fill_lim_pct <- max(abs(scatter_df$pct_clin), na.rm = TRUE) * c(-1, 1)
  pct_breaks   <- waiver()
  pct_labels   <- waiver()
} else {
  fill_lim_pct <- c(-pct_scatter_lim, pct_scatter_lim)
  lo_lab <- if (pct_scatter_lim >= 100) "-100%" else sprintf("<-%g%%", pct_scatter_lim)
  hi_lab <- sprintf(">%g%%", pct_scatter_lim)
  pct_breaks <- c(-pct_scatter_lim, -pct_scatter_lim / 2, 0,
                   pct_scatter_lim / 2,  pct_scatter_lim)
  pct_labels <- c(lo_lab,
                  sprintf("-%g%%", pct_scatter_lim / 2),
                  "0%",
                  sprintf("+%g%%", pct_scatter_lim / 2),
                  hi_lab)
}

p_scatter <- ggplot(scatter_df, aes(res_fut, eir, fill = delta_clin)) +
  geom_point(shape = 21, size = 6, colour = "grey25", stroke = 0.5) +
  ggrepel::geom_text_repel(aes(label = region), size = 3, colour = "grey20",
                           box.padding = 0.5, seed = 1, max.overlaps = Inf) +
  scale_y_log10() +
  cmocean::scale_fill_cmocean(name = "balance", limits = fill_lim) +
  theme_minimal(base_size = 11) +
  labs(
    x = "Future pyrethroid resistance (mean over distribution years)",
    y = "Baseline EIR (log scale)",
    fill = "ATN - Pyr-CFP\nclinical cases\n(/1,000 / yr)",
    title = paste0(country_name, " — ATN vs Pyr-CFP by resistance and transmission intensity"),
    caption = paste0(
      "Each point = one admin-1 region. Fill = annual clinical incidence under ATN minus Pyr-CFP",
      " (per 1,000 / yr over the ", meta$n_future_years, "-yr window).\n",
      "Blue = ATN averts more than Pyr-CFP (better); red = ATN averts less (worse).",
      " x = arithmetic mean future resistance (", min(fut_years), "-", max(fut_years), ")."))

p_scatter_pct <- ggplot(scatter_df, aes(res_fut, eir, fill = pct_clin)) +
  geom_point(shape = 21, size = 6, colour = "grey25", stroke = 0.5) +
  ggrepel::geom_text_repel(aes(label = region), size = 3, colour = "grey20",
                           box.padding = 0.5, seed = 1, max.overlaps = Inf) +
  scale_y_log10() +
  cmocean::scale_fill_cmocean(name = "balance", limits = fill_lim_pct,
                              oob = scales::squish,
                              breaks = pct_breaks, labels = pct_labels) +
  theme_minimal(base_size = 11) +
  labs(
    x = "Future pyrethroid resistance (mean over distribution years)",
    y = "Baseline EIR (log scale)",
    fill = "ATN vs Pyr-CFP\nclinical cases\n(% change)",
    title = paste0(country_name, " — ATN vs Pyr-CFP by resistance and transmission intensity (% change)"),
    caption = paste0(
      "Each point = one admin-1 region. Fill = (ATN - Pyr-CFP) / Pyr-CFP x 100%",
      " over the ", meta$n_future_years, "-yr window.\n",
      "Blue = ATN has fewer cases than Pyr-CFP (better); red = ATN has more cases (worse).",
      " x = arithmetic mean future resistance (", min(fut_years), "-", max(fut_years), ")."))

out_scatter     <- sprintf("dev/outputs/%s_resistance_eir_scatter.png",     tolower(COUNTRY_ISO))
out_scatter_pct <- sprintf("dev/outputs/%s_resistance_eir_scatter_pct.png", tolower(COUNTRY_ISO))
ggsave(out_scatter,     p_scatter,     width = 9, height = 7, dpi = 150)
ggsave(out_scatter_pct, p_scatter_pct, width = 9, height = 7, dpi = 150)
message(sprintf("Saved: %s, %s", out_scatter, out_scatter_pct))

# ---------------------------------------------------------------------
# 3. Past-window consistency check
# ---------------------------------------------------------------------
past_check <- df |>
  filter(year_rel < 0) |>
  group_by(region, timestep) |>
  summarise(n_distinct_pfpr = n_distinct(round(pfpr2to10, 6)), .groups = "drop")
if (any(past_check$n_distinct_pfpr > 1)) {
  warning("Arms diverge in the past window (expected: IBM stochasticity + Erlang vs discrete-delay ODE).")
} else {
  message("Verification OK: all arms identical in the past window.")
}

print(p_clin_series)
print(p_prev)
print(p_maps)
