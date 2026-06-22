# =====================================================================
# c24med_projection_plots.R
#   Reads dev/outputs/<iso>_c24med_projection_results.rds and produces
#   prevalence/incidence geofacet timeseries + cases-averted choropleths
#   for all 7 arms (none / pyr / pyr_pbo / pyr_cfp / atn / pyr_atn /
#   pyr_cfp_atn).
#   Change COUNTRY_ISO to plot any available c24med run.
# =====================================================================

library(dplyr); library(tidyr); library(ggplot2); library(zoo)
library(sf); library(patchwork)

# ---------------------------------------------------------------------
# 0. Config
# ---------------------------------------------------------------------
COUNTRY_ISO  <- "MLI"
RESULTS_FILE <- sprintf("dev/outputs/%s_c24med_projection_results.rds", tolower(COUNTRY_ISO))
SITE_FILE    <- sprintf("dev/site_files/without_split/%s.rds", COUNTRY_ISO)

obj  <- readRDS(RESULTS_FILE)
df   <- obj$results
shape <- obj$shape
meta <- obj$meta

country_name <- if (!is.null(meta$country_name)) meta$country_name else COUNTRY_ISO
SHAPE_KEY    <- "name_1"

# ---------------------------------------------------------------------
# 1. Colour assignments — Paul Tol "muted" qualitative palette
#
#  Full tol-muted palette (9 colours + 1 pale "bad-data" grey):
#    rose:   #CC6677
#    indigo: #332288
#    sand:   #DDCC77
#    green:  #117733
#    cyan:   #88CCEE
#    wine:   #882255
#    teal:   #44AA99
#    olive:  #999933
#    purple: #AA4499
#    pale:   #DDDDDD  (use for missing / "bad" data only)
#
#  Arm assignments (adjust hex codes here to remap colours):
# ---------------------------------------------------------------------
arm_labels <- c(
  none        = "No future nets",
  pyr         = "Future Pyr",
  pyr_pbo     = "Future Pyr-PBO",
  pyr_cfp     = "Future Pyr-CFP",
  atn         = "Future ATN",
  pyr_atn     = "Future Pyr-ATN",
  pyr_cfp_atn = "Future Pyr-CFP-ATN"
)

arm_cols <- c(
  "No future nets"      = "#DDDDDD",  # pale   — baseline / no intervention
  "Future Pyr"          = "#88CCEE",  # cyan   — pyrethroid only
  "Future Pyr-PBO"      = "#44AA99",  # teal   — pyrethroid + PBO
  "Future Pyr-CFP"      = "#117733",  # green  — pyrethroid + CFP
  "Future ATN"          = "#DDCC77",  # sand   — antimalarial net (no insecticide)
  "Future Pyr-ATN"      = "#CC6677",  # rose   — pyrethroid + antimalarial
  "Future Pyr-CFP-ATN"  = "#882255"   # wine   — pyrethroid-CFP + antimalarial
)

# ---------------------------------------------------------------------
# 2. Prep windows
# ---------------------------------------------------------------------
df_win <- df |>
  filter(year_rel >= 0, year_rel <= meta$n_future_years) |>
  mutate(arm_f = factor(arm_labels[arm], levels = arm_labels))

missing_regions <- setdiff(shape[[SHAPE_KEY]], unique(df$region))
excl_note <- if (length(missing_regions) > 0L) {
  paste0(" Excluded (no successful run): ", paste(missing_regions, collapse = ", "), ".")
} else {
  ""
}

# ---------------------------------------------------------------------
# 3. Geofacet timeseries (prevalence + clinical incidence)
# ---------------------------------------------------------------------
country_grid_raw <- geofacet::grid_auto(shape, names = SHAPE_KEY, seed = 1)
names(country_grid_raw)[names(country_grid_raw) == paste0("name_", SHAPE_KEY)] <- "name"
country_grid_raw$code <- country_grid_raw$name
present_regions <- unique(df_win$region)
country_grid <- country_grid_raw[country_grid_raw$code %in% present_regions,
                                 c("row", "col", "code", "name")]

df_plot <- df |>
  filter(year_rel >= -3) |>
  mutate(
    arm_f     = factor(arm_labels[arm], levels = arm_labels),
    cal_year  = year_rel + meta$future_yr0,
    clin_rate = clin_inc / meta$human_pop * 1000
  ) |>
  group_by(region, arm_f) |>
  arrange(cal_year, .by_group = TRUE) |>
  mutate(
    pfpr_roll      = zoo::rollmean(pfpr2to10 * 100, k = 365L, fill = NA, align = "center"),
    clin_rate_roll = zoo::rollmean(clin_rate,        k = 365L, fill = NA, align = "center")
  ) |>
  ungroup()

vline_df     <- data.frame(xintercept = meta$future_yr0)
plot_caption <- paste0(
  "Dashed line = future distribution start (", meta$future_yr0, ").",
  " Thick line = 365-day rolling mean.\n",
  "'No future nets': historical ITNs only, decaying from ", meta$future_yr0,
  " (no new distributions).", excl_note,
  "\nNet efficacy: Churcher 2024 median estimates."
)

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
       title   = paste0(country_name,
                        " — projected prevalence by region and future net scenario (c24med)"),
       caption = plot_caption)

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
       title   = paste0(country_name,
                        " — projected clinical incidence by region and future net scenario (c24med)"),
       caption = plot_caption)

out_prev <- sprintf("dev/outputs/%s_c24med_prevalence_facet.png",  tolower(COUNTRY_ISO))
out_clin <- sprintf("dev/outputs/%s_c24med_incidence_facet.png",   tolower(COUNTRY_ISO))
ggsave(out_prev, p_prev,         width = 13, height = 9, dpi = 150)
ggsave(out_clin, p_clin_series,  width = 13, height = 9, dpi = 150)
message(sprintf("Saved: %s, %s", out_prev, out_clin))

# ---------------------------------------------------------------------
# 4. Clinical cases averted vs no-nets (future window, per 1,000 pop / yr)
#    Six intervention arms shown as two rows of three maps.
# ---------------------------------------------------------------------
per1000_yr <- function(total) total / meta$human_pop * 1000 / meta$n_future_years

tot <- df_win |>
  group_by(region, arm) |>
  summarise(cases = sum(clin_inc), .groups = "drop") |>
  pivot_wider(names_from = arm, values_from = cases) |>
  mutate(
    pyr_averted         = per1000_yr(none - pyr),
    pyr_pbo_averted     = per1000_yr(none - pyr_pbo),
    pyr_cfp_averted     = per1000_yr(none - pyr_cfp),
    atn_averted         = per1000_yr(none - atn),
    pyr_atn_averted     = per1000_yr(none - pyr_atn),
    pyr_cfp_atn_averted = per1000_yr(none - pyr_cfp_atn)
  )

map_df <- shape |> left_join(tot, by = setNames("region", SHAPE_KEY))

# Shared fill limits across all six averted maps so colours are comparable.
averted_cols <- c("pyr_averted", "pyr_pbo_averted", "pyr_cfp_averted",
                  "atn_averted", "pyr_atn_averted", "pyr_cfp_atn_averted")
averted_rng  <- range(unlist(tot[averted_cols]), na.rm = TRUE)

make_averted_map <- function(fill_col, title) {
  ggplot(map_df) +
    geom_sf(aes(fill = .data[[fill_col]]), colour = "white", linewidth = 0.2) +
    scale_fill_viridis_c(option = "D", limits = averted_rng) +
    theme_void(base_size = 10) +
    labs(title = title, fill = "Averted /\n1,000 / yr")
}

# Add a Bamako bounding-box rectangle to the main map and overlay a zoomed
# inset in the north-west (top-left) corner of the panel.
# Adjust `expand` (degrees) and inset corner coords (0-1 npc) to taste.
add_bamako_inset <- function(p_main, fill_col, fill_lims) {
  bamako_sf <- map_df[map_df[[SHAPE_KEY]] == "Bamako", ]

  # Expand the bounding box slightly so the rectangle is visible on the main map.
  bb <- sf::st_bbox(bamako_sf)
  expand <- 0.5   # degrees; increase if the highlight box looks too tight
  bb["xmin"] <- bb["xmin"] - expand
  bb["ymin"] <- bb["ymin"] - expand
  bb["xmax"] <- bb["xmax"] + expand
  bb["ymax"] <- bb["ymax"] + expand
  bb_rect <- sf::st_as_sfc(bb)

  p_rect <- p_main +
    geom_sf(data = bb_rect, fill = NA, colour = "black",
            linewidth = 0.7, inherit.aes = FALSE)

  # Inset: Bamako polygon with matching fill scale, black border, no legend.
  p_inset <- ggplot(bamako_sf) +
    geom_sf(aes(fill = .data[[fill_col]]), colour = "grey30", linewidth = 0.5) +
    scale_fill_viridis_c(option = "D", limits = fill_lims) +
    theme_void(base_size = 8) +
    theme(
      legend.position = "none",
      panel.border    = element_rect(colour = "black", fill = NA, linewidth = 0.8)
    )

  # Overlay in the NW corner; coords are 0-1 relative to the map panel.
  # Adjust left/bottom/right/top if the inset needs repositioning.
  p_rect + patchwork::inset_element(
    p_inset,
    left = 0, bottom = 0.60, right = 0.35, top = 1.0,
    align_to = "panel"
  )
}

m_pyr         <- add_bamako_inset(make_averted_map("pyr_averted",         "Pyr"),
                                  "pyr_averted",         averted_rng)
m_pyr_pbo     <- add_bamako_inset(make_averted_map("pyr_pbo_averted",     "Pyr-PBO"),
                                  "pyr_pbo_averted",     averted_rng)
m_pyr_cfp     <- add_bamako_inset(make_averted_map("pyr_cfp_averted",     "Pyr-CFP"),
                                  "pyr_cfp_averted",     averted_rng)
m_atn         <- add_bamako_inset(make_averted_map("atn_averted",         "ATN"),
                                  "atn_averted",         averted_rng)
m_pyr_atn     <- add_bamako_inset(make_averted_map("pyr_atn_averted",     "Pyr-ATN"),
                                  "pyr_atn_averted",     averted_rng)
m_pyr_cfp_atn <- add_bamako_inset(make_averted_map("pyr_cfp_atn_averted", "Pyr-CFP-ATN"),
                                  "pyr_cfp_atn_averted", averted_rng)

maps_caption <- paste0(
  "'No future nets' = historical ITNs decaying from ", meta$future_yr0,
  " with no replacements.\n",
  "All maps share the same fill scale. Net efficacy: Churcher 2024 median estimates.",
  excl_note
)

p_maps <- (m_pyr | m_pyr_pbo | m_pyr_cfp) / (m_atn | m_pyr_atn | m_pyr_cfp_atn) +
  patchwork::plot_annotation(
    title   = paste0(country_name,
                     " — clinical cases averted vs no-nets over ", meta$n_future_years,
                     "-year projection (per 1,000 pop / yr, c24med)"),
    caption = maps_caption
  )

out_maps <- sprintf("dev/outputs/%s_c24med_averted_maps.png", tolower(COUNTRY_ISO))
ggsave(out_maps, p_maps, width = 12, height = 8, dpi = 150)
message(sprintf("Saved: %s", out_maps))

# ---------------------------------------------------------------------
# 5. Past-window consistency check
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

print(p_prev)
print(p_clin_series)
print(p_maps)
