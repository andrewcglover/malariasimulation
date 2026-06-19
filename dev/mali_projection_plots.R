# =====================================================================
# mali_projection_plots.R
#   Reads dev/outputs/mali_projection_results.rds and produces:
#     (1) prevalence-over-time, geographically facetted by region (geofacet)
#     (2) three clinical-cases-averted choropleths:
#           Pyr-CFP vs none, ATN vs none, ATN additional vs Pyr-CFP
#   Run from the fork root after mali_projection_run.R has completed.
# =====================================================================

library(dplyr); library(tidyr); library(ggplot2)
library(sf); library(patchwork)

obj  <- readRDS("dev/outputs/mali_projection_results.rds")
df   <- obj$results
shape <- obj$shape    # sf object for MLI admin-1 (level_1 polygons)
meta <- obj$meta

# CONFIRMED: shape column holding region names is "name_1"
SHAPE_KEY <- "name_1"

arm_labels <- c(none    = "No future nets",
                cfp     = "Future Pyr-CFP",
                atn     = "Future ATN",
                pyr_atn = "Future Pyr-ATN")
arm_cols   <- c("No future nets" = "grey50",
                "Future Pyr-CFP" = "#009988",
                "Future ATN"     = "#EE7733",
                "Future Pyr-ATN" = "#CC3311")

# ---------------------------------------------------------------------
# 1. Geographic facets of prevalence over the future window
# ---------------------------------------------------------------------
df_win <- df |>
  filter(year_rel >= 0, year_rel <= meta$n_future_years) |>
  mutate(arm_f = factor(arm_labels[arm], levels = arm_labels))

# Data-driven geofacet grid from the shape centroids.
# grid_auto() needs a df with code, name, lon, lat columns.
cent     <- suppressWarnings(sf::st_centroid(shape))
xy       <- as.data.frame(sf::st_coordinates(cent))
auto_grid <- data.frame(
  name = shape[[SHAPE_KEY]],
  code = shape[[SHAPE_KEY]],
  lon  = xy$X,
  lat  = xy$Y
)

# Try grid_auto; fall back to the hand-tuned grid if unavailable or layout looks wrong.
mali_grid <- tryCatch(
  geofacet::grid_auto(auto_grid, codes = "code", names = "name",
                      seed = 42),
  error = function(e) NULL
)

if (is.null(mali_grid)) {
  # Hand-tuned fallback — north at top. Region names MUST match unique(df$region)
  # exactly (i.e. the site file spellings: Timbuktu, Ségou, etc.).
  mali_grid <- data.frame(
    row  = c(1,    1,    1,          2,     2,         3,         3,       4,        4),
    col  = c(2,    3,    5,          3,     4,         1,         2,       2,        3),
    code = c("Timbuktu","Gao","Kidal","Mopti","Ségou","Kayes","Koulikoro","Bamako","Sikasso"),
    name = c("Timbuktu","Gao","Kidal","Mopti","Ségou","Kayes","Koulikoro","Bamako","Sikasso")
  )
}

p_prev <- ggplot(df_win, aes(year_rel, pfpr2to10 * 100, colour = arm_f)) +
  geom_line(linewidth = 0.6) +
  scale_colour_manual(values = arm_cols) +
  geofacet::facet_geo(~ region, grid = mali_grid) +
  theme_minimal(base_size = 11) +
  labs(x = "Years post first future distribution",
       y = expression(italic(Pf) * PR[2 - 10] * " (%)"),
       colour = "",
       title = "Mali — projected prevalence by region and future net scenario")

ggsave("dev/outputs/mali_prevalence_geofacet.png", p_prev,
       width = 11, height = 8, dpi = 150)

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

# Join to admin-1 sf polygons on the shared name_1 / region key.
# setNames("region", SHAPE_KEY) -> c(name_1 = "region"):
#   left_join(shape, tot, by = c(name_1 = "region"))
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
m_cfp <- make_map("cfp_averted",    "Pyr-CFP vs none") +
           scale_fill_viridis_c(limits = rng, option = "D")
m_atn <- make_map("atn_averted",    "ATN vs none") +
           scale_fill_viridis_c(limits = rng, option = "D")
m_add <- make_map("atn_additional", "ATN additional (vs Pyr-CFP)", diverging = TRUE)

p_maps <- (m_cfp | m_atn | m_add) +
  patchwork::plot_annotation(
    title = "Mali — clinical cases averted over 6-year projection (per 1,000 pop / yr)"
  )

ggsave("dev/outputs/mali_averted_maps.png", p_maps,
       width = 14, height = 5, dpi = 150)

message("Saved: dev/outputs/mali_prevalence_geofacet.png, dev/outputs/mali_averted_maps.png")

# ---------------------------------------------------------------------
# 3. Past-window consistency check
#
# Two expected sources of arm divergence even before future_start_day:
#   (a) IBM stochasticity: each arm is an independent run_simulation() call
#       with no shared seed, so the human IBM drifts independently.
#   (b) ODE formulation: atn/pyr_atn arms use the Erlang EIP chain
#       (spor_len stages) while none/cfp use the discrete-delay model;
#       these converge to tolerance but are not bit-identical (see CLAUDE.md §5).
# A warning here reflects structural design, not a parameter bug.
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
print(p_maps)
