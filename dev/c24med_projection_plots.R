# =====================================================================
# c24med_projection_plots.R
#   Reads dev/outputs/<iso>_c24med_projection_results_hl<tag>.rds and produces:
#     (1) prevalence / incidence geofacet timeseries
#     (2) cases-averted vs no-nets bar chart (geofacet)
#     (3) cases-averted choropleths (6 arms, 2x3 layout)
#     (4) % ATN-exposed mosquitoes geofacet (ATN-bearing arms only)
#   Change COUNTRY_ISO / HL_TAG to plot any available c24med run.
# =====================================================================

library(dplyr); library(tidyr); library(ggplot2); library(zoo)
library(sf); library(patchwork)

# ---------------------------------------------------------------------
# 0. Config
# ---------------------------------------------------------------------
COUNTRY_ISO  <- "MLI"
HL_TAG       <- "hl2p64"
RESULTS_FILE <- sprintf("dev/outputs/%s_c24med_projection_results_%s.rds",
                        tolower(COUNTRY_ISO), HL_TAG)
SITE_FILE    <- sprintf("dev/site_files/without_split/%s.rds", COUNTRY_ISO)

obj   <- readRDS(RESULTS_FILE)
df    <- obj$results
shape <- obj$shape
meta  <- obj$meta

SHAPE_KEY <- "name_1"

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
  none        = "No nets",
  pyr         = "Pyr",
  pyr_pbo     = "Pyr-PBO",
  pyr_cfp     = "Pyr-CFP",
  atn         = "ATN",
  pyr_atn     = "Pyr-ATN",
  pyr_cfp_atn = "Pyr-CFP-ATN"
)

arm_cols <- c(
  "No nets"       = "#DDDDDD",  # pale   — no intervention baseline
  "Pyr"           = "#882255",  # wine
  "Pyr-PBO"       = "#CC6677",  # rose
  "Pyr-CFP"       = "#DDCC77",  # sand
  "ATN"           = "#117733",  # green
  "Pyr-ATN"       = "#44AA99",  # teal
  "Pyr-CFP-ATN"   = "#88CCEE"   # cyan
)

# ---------------------------------------------------------------------
# 2. Prep windows & grid
# ---------------------------------------------------------------------
df_win <- df |>
  filter(year_rel >= 0, year_rel <= meta$n_future_years) |>
  mutate(arm_f = factor(arm_labels[arm], levels = arm_labels))

present_regions  <- unique(df_win$region)
missing_regions  <- setdiff(shape[[SHAPE_KEY]], present_regions)

country_grid_raw <- geofacet::grid_auto(shape, names = SHAPE_KEY, seed = 1)
names(country_grid_raw)[names(country_grid_raw) == paste0("name_", SHAPE_KEY)] <- "name"
country_grid_raw$code <- country_grid_raw$name
country_grid <- country_grid_raw[country_grid_raw$code %in% present_regions,
                                 c("row", "col", "code", "name")]

vline_df <- data.frame(xintercept = meta$future_yr0)

# output filename helper
out_file <- function(suffix) {
  sprintf("dev/outputs/%s_c24med_%s_%s.png", tolower(COUNTRY_ISO), suffix, HL_TAG)
}

# ---------------------------------------------------------------------
# 3. Geofacet timeseries — prevalence & clinical incidence
# ---------------------------------------------------------------------
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

p_prev <- ggplot(df_plot, aes(cal_year, pfpr2to10 * 100, colour = arm_f)) +
  geom_vline(data = vline_df, aes(xintercept = xintercept),
             linetype = "dashed", colour = "grey40", linewidth = 0.4) +
  geom_line(linewidth = 0.4, alpha = 0.2) +
  geom_line(aes(y = pfpr_roll), linewidth = 0.8, na.rm = TRUE) +
  scale_colour_manual(values = arm_cols) +
  geofacet::facet_geo(~ region, grid = country_grid) +
  theme_minimal(base_size = 15) +
  labs(x = "Year",
       y = expression(italic(Pf) * PR[2 - 10] * " (%)"),
       colour = "")

p_clin_series <- ggplot(df_plot, aes(cal_year, clin_rate, colour = arm_f)) +
  geom_vline(data = vline_df, aes(xintercept = xintercept),
             linetype = "dashed", colour = "grey40", linewidth = 0.4) +
  geom_line(linewidth = 0.4, alpha = 0.2) +
  geom_line(aes(y = clin_rate_roll), linewidth = 0.8, na.rm = TRUE) +
  scale_colour_manual(values = arm_cols) +
  geofacet::facet_geo(~ region, grid = country_grid) +
  theme_minimal(base_size = 15) +
  labs(x = "Year",
       y = "Clinical incidence (per 1,000 / day, all ages)",
       colour = "")

ggsave(out_file("prevalence_facet"), p_prev,        width = 13, height = 9, dpi = 300)
ggsave(out_file("incidence_facet"),  p_clin_series, width = 13, height = 9, dpi = 300)
message(sprintf("Saved: %s, %s", out_file("prevalence_facet"), out_file("incidence_facet")))

# ---------------------------------------------------------------------
# 4. Bar chart — cases averted vs no-nets (geofacet)
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

averted_long <- tot |>
  select(region, ends_with("_averted")) |>
  pivot_longer(-region, names_to = "arm_col", values_to = "averted") |>
  mutate(
    arm   = sub("_averted$", "", arm_col),
    arm_f = factor(arm_labels[arm], levels = arm_labels)
  )

p_bar <- ggplot(averted_long, aes(x = arm_f, y = averted, fill = arm_f)) +
  geom_col(width = 0.75) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey40") +
  scale_fill_manual(values = arm_cols) +
  geofacet::facet_geo(~ region, grid = country_grid, scales = "fixed") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x     = element_blank(),
    axis.ticks.x    = element_blank(),
    legend.position = "bottom"
  ) +
  guides(fill = guide_legend(nrow = 1)) +
  labs(x = "", y = "Cases averted vs no nets (per 1,000 / yr)", fill = "")

ggsave(out_file("averted_bar"), p_bar, width = 13, height = 9, dpi = 300)
message(sprintf("Saved: %s", out_file("averted_bar")))

# ---------------------------------------------------------------------
# 5. Averted choropleths (2 x 3), with Bamako inset for MLI only
# ---------------------------------------------------------------------
map_df <- shape |> left_join(tot, by = setNames("region", SHAPE_KEY))

averted_cols <- c("pyr_averted", "pyr_pbo_averted", "pyr_cfp_averted",
                  "atn_averted", "pyr_atn_averted", "pyr_cfp_atn_averted")
averted_rng  <- range(unlist(tot[averted_cols]), na.rm = TRUE)

make_averted_map <- function(fill_col, title) {
  ggplot(map_df) +
    geom_sf(aes(fill = .data[[fill_col]]), colour = "white", linewidth = 0.2) +
    scale_fill_viridis_c(option = "D", limits = averted_rng) +
    theme_void(base_size = 14) +
    labs(title = title, fill = "Averted /\n1,000 / yr")
}

# Bamako inset — only applied for MLI (where Bamako is the small capital district).
# For other countries, make_averted_map() output is used directly.
# Padding is 30% of Bamako's own width/height so the rectangle on the main map
# and the inset zoom window have the same proportional margin around the district.
# Inset panel position: left/bottom/right/top in 0-1 npc relative to map panel.
add_bamako_inset <- function(p_main, fill_col, fill_lims) {
  bamako_sf <- map_df[map_df[[SHAPE_KEY]] == "Bamako", ]
  bb        <- sf::st_bbox(bamako_sf)

  pad_frac <- 0.30   # 30% of Bamako's own extent on each side
  expand_x <- as.numeric(bb["xmax"] - bb["xmin"]) * pad_frac
  expand_y <- as.numeric(bb["ymax"] - bb["ymin"]) * pad_frac

  xlim <- c(bb["xmin"] - expand_x, bb["xmax"] + expand_x)
  ylim <- c(bb["ymin"] - expand_y, bb["ymax"] + expand_y)

  bb["xmin"] <- xlim[1]; bb["xmax"] <- xlim[2]
  bb["ymin"] <- ylim[1]; bb["ymax"] <- ylim[2]
  bb_rect <- sf::st_as_sfc(bb)

  p_rect <- p_main +
    geom_sf(data = bb_rect, fill = NA, colour = "black",
            linewidth = 0.7, inherit.aes = FALSE)

  p_inset <- ggplot(bamako_sf) +
    geom_sf(aes(fill = .data[[fill_col]]), colour = "grey30", linewidth = 0.5) +
    scale_fill_viridis_c(option = "D", limits = fill_lims) +
    coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
    theme_void(base_size = 11) +
    theme(
      legend.position = "none",
      panel.border    = element_rect(colour = "black", fill = NA, linewidth = 0.8)
    )

  p_rect + patchwork::inset_element(
    p_inset,
    left = 0, bottom = 0.72, right = 0.26, top = 1.0,
    align_to = "panel"
  )
}

wrap_map <- function(p_main, fill_col, fill_lims) {
  if (COUNTRY_ISO == "MLI") {
    add_bamako_inset(p_main, fill_col, fill_lims)
  } else {
    p_main
  }
}

m_pyr         <- wrap_map(make_averted_map("pyr_averted",         "Pyr"),         "pyr_averted",         averted_rng)
m_pyr_pbo     <- wrap_map(make_averted_map("pyr_pbo_averted",     "Pyr-PBO"),     "pyr_pbo_averted",     averted_rng)
m_pyr_cfp     <- wrap_map(make_averted_map("pyr_cfp_averted",     "Pyr-CFP"),     "pyr_cfp_averted",     averted_rng)
m_atn         <- wrap_map(make_averted_map("atn_averted",         "ATN"),         "atn_averted",         averted_rng)
m_pyr_atn     <- wrap_map(make_averted_map("pyr_atn_averted",     "Pyr-ATN"),     "pyr_atn_averted",     averted_rng)
m_pyr_cfp_atn <- wrap_map(make_averted_map("pyr_cfp_atn_averted", "Pyr-CFP-ATN"), "pyr_cfp_atn_averted", averted_rng)

p_maps <- (m_pyr | m_pyr_pbo | m_pyr_cfp) / (m_atn | m_pyr_atn | m_pyr_cfp_atn) +
  patchwork::plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(out_file("averted_maps"), p_maps, width = 12, height = 9, dpi = 300)
message(sprintf("Saved: %s", out_file("averted_maps")))

# ---------------------------------------------------------------------
# 6. % ATN-exposed mosquitoes over time (geofacet, ATN-bearing arms only)
#    Summed across all species and all adult compartments (Sv + Ev + Iv).
# ---------------------------------------------------------------------
exp_cols <- grep("^(Sv|Ev|Iv)_exposed_",             names(df), value = TRUE)
all_cols <- grep("^(Sv|Ev|Iv)_(exposed|unexposed)_", names(df), value = TRUE)

atn_arms <- c("atn", "pyr_atn", "pyr_cfp_atn")

df_exp <- df |>
  filter(arm %in% atn_arms, year_rel >= -3) |>
  mutate(
    arm_f      = factor(arm_labels[arm], levels = arm_labels),
    cal_year   = year_rel + meta$future_yr0,
    mosq_exp   = rowSums(across(all_of(exp_cols))),
    mosq_total = rowSums(across(all_of(all_cols))),
    pct_exp    = if_else(mosq_total > 0, mosq_exp / mosq_total * 100, NA_real_)
  ) |>
  group_by(region, arm_f) |>
  arrange(cal_year, .by_group = TRUE) |>
  mutate(pct_exp_roll = zoo::rollmean(pct_exp, k = 365L, fill = NA, align = "center")) |>
  ungroup()

p_pct_exp <- ggplot(df_exp, aes(cal_year, pct_exp, colour = arm_f)) +
  geom_vline(data = vline_df, aes(xintercept = xintercept),
             linetype = "dashed", colour = "grey40", linewidth = 0.4) +
  geom_line(linewidth = 0.4, alpha = 0.2) +
  geom_line(aes(y = pct_exp_roll), linewidth = 0.8, na.rm = TRUE) +
  scale_colour_manual(values = arm_cols) +
  geofacet::facet_geo(~ region, grid = country_grid) +
  theme_minimal(base_size = 15) +
  labs(x = "Year",
       y = "ATN-exposed adult mosquitoes (%)",
       colour = "")

ggsave(out_file("pct_atn_exposed"), p_pct_exp, width = 13, height = 9, dpi = 300)
message(sprintf("Saved: %s", out_file("pct_atn_exposed")))

# ---------------------------------------------------------------------
# 7. Past-window consistency check
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
print(p_bar)
print(p_maps)
print(p_pct_exp)
