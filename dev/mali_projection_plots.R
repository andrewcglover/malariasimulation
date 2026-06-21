# =====================================================================
# mali_projection_plots.R
#   Reads dev/outputs/mali_projection_results.rds and produces:
#     (1) prevalence-over-time, geographically facetted by region (geofacet)
#     (2) three clinical-cases-averted choropleths:
#           Pyr-CFP vs none, ATN vs none, ATN additional vs Pyr-CFP
#   Run from the fork root after mali_projection_run.R has completed.
# =====================================================================

library(dplyr); library(tidyr); library(ggplot2); library(zoo)
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

# Build geofacet grid via grid_auto() — uses geogrid's hexagonal-regular
# layout algorithm so the resulting row/col positions are guaranteed
# compatible with facet_geo's internal matching.  Takes names from the
# shape object (encoding-safe), then filters to regions present in the data
# (drops Bamako which failed the ODE solver).
mali_grid_raw <- geofacet::grid_auto(shape, names = SHAPE_KEY, seed = 1)
# grid_auto prefixes "name_" to the column name; rename to plain "name" / add "code"
names(mali_grid_raw)[names(mali_grid_raw) == paste0("name_", SHAPE_KEY)] <- "name"
mali_grid_raw$code <- mali_grid_raw$name
present_regions <- unique(df_win$region)
mali_grid <- mali_grid_raw[mali_grid_raw$code %in% present_regions,
                            c("row", "col", "code", "name")]

# df_plot: 3 years history + future window, with calendar year and rolling means.
# Kept separate from df_win so the cases-averted sum (section 2) is unaffected.
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
    pfpr_roll     = zoo::rollmean(pfpr2to10 * 100, k = 365L, fill = NA, align = "center"),
    clin_rate_roll = zoo::rollmean(clin_rate,       k = 365L, fill = NA, align = "center")
  ) |>
  ungroup()

vline_df <- data.frame(xintercept = meta$future_yr0)
plot_caption <- paste0("Dashed line = future distribution start (", meta$future_yr0, ").",
                       " Thick line = 365-day rolling mean.",
                       "\n'No future nets': historical ITNs only, decaying from ", meta$future_yr0,
                       " (no new distributions). Bamako excluded: ODE solver failure.")

p_clin_series <- ggplot(df_plot, aes(cal_year, clin_rate, colour = arm_f)) +
  geom_vline(data = vline_df, aes(xintercept = xintercept),
             linetype = "dashed", colour = "grey40", linewidth = 0.4) +
  geom_line(linewidth = 0.4, alpha = 0.2) +
  geom_line(aes(y = clin_rate_roll), linewidth = 0.8, na.rm = TRUE) +
  scale_colour_manual(values = arm_cols) +
  geofacet::facet_geo(~ region, grid = mali_grid) +
  theme_minimal(base_size = 11) +
  labs(x = "Year",
       y = "Clinical incidence (per 1,000 / day, all ages)",
       colour = "",
       title = "Mali — projected clinical incidence by region and future net scenario",
       caption = plot_caption)

p_prev <- ggplot(df_plot, aes(cal_year, pfpr2to10 * 100, colour = arm_f)) +
  geom_vline(data = vline_df, aes(xintercept = xintercept),
             linetype = "dashed", colour = "grey40", linewidth = 0.4) +
  geom_line(linewidth = 0.4, alpha = 0.2) +
  geom_line(aes(y = pfpr_roll), linewidth = 0.8, na.rm = TRUE) +
  scale_colour_manual(values = arm_cols) +
  geofacet::facet_geo(~ region, grid = mali_grid) +
  theme_minimal(base_size = 11) +
  labs(x = "Year",
       y = expression(italic(Pf) * PR[2 - 10] * " (%)"),
       colour = "",
       title = "Mali — projected prevalence by region and future net scenario",
       caption = plot_caption)

ggsave("dev/outputs/mali_prevalence_facet.png",  p_prev,       width = 11, height = 8, dpi = 150)
ggsave("dev/outputs/mali_incidence_facet.png",   p_clin_series, width = 11, height = 8, dpi = 150)

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
m_cfp <- make_map("cfp_averted",    "Pyr-CFP vs\nno future nets") +
           scale_fill_viridis_c(limits = rng, option = "D")
m_atn <- make_map("atn_averted",    "ATN vs\nno future nets") +
           scale_fill_viridis_c(limits = rng, option = "D")
m_add <- make_map("atn_additional", "ATN additional\n(vs Pyr-CFP)", diverging = TRUE)

p_maps <- (m_cfp | m_atn | m_add) +
  patchwork::plot_annotation(
    title = "Mali — clinical cases averted over 6-year projection (per 1,000 pop / yr)",
    caption = paste0("'No future nets' = historical ITNs decaying from 2025 with no replacements;",
                     " incidence rebounds toward unprotected equilibrium.\n",
                     "High averted counts reflect this rebound, not the current protected baseline.",
                     " Bamako excluded (ODE solver failure).")
  )

ggsave("dev/outputs/mali_averted_maps.png", p_maps,
       width = 14, height = 5, dpi = 150)

message("Saved: dev/outputs/mali_prevalence_facet.png, dev/outputs/mali_incidence_facet.png, dev/outputs/mali_averted_maps.png")

# ---------------------------------------------------------------------
# 2b. Resistance x baseline-EIR scatter, coloured by ATN - Pyr-CFP
#     change in annual clinical incidence (per 1,000 / yr).
#
#   x = future pyrethroid resistance, ARITHMETIC mean over the future
#       distribution years (future_yr0 .. future_yr0 + n_future_years - 1).
#       The model itself applies per-distribution-year resistance
#       (build_future_schedule); this mean is only a per-region summary for
#       the axis. Resistance is a bounded survival proportion, not a
#       multiplicative rate, so arithmetic (not geometric) is the natural
#       summary; over this 0.84-0.96 range the two differ by < 1e-3 anyway.
#   y = baseline (calibration) EIR from the site file, log10 scale
#       (spans ~6 to ~500 across regions).
#   fill = ATN cases - Pyr-CFP cases, per 1,000 / yr (cmocean "balance",
#       symmetric about 0):  RED  = ATN worse (more cases),
#                            BLUE = ATN better (fewer cases).
# ---------------------------------------------------------------------
site_obj <- readRDS("dev/site_files/without_split/MLI.rds")

fut_years <- meta$future_yr0 + 0:(meta$n_future_years - 1L)
res_fut <- site_obj$vectors$pyrethroid_resistance |>
  filter(year %in% fut_years) |>
  group_by(name_1) |>
  summarise(res_fut = mean(pyrethroid_resistance), .groups = "drop")

eir_base <- site_obj$eir |>
  filter(sp == "pf") |>
  select(name_1, eir)

# ATN - Pyr-CFP change in clinical incidence (positive = ATN has MORE cases).
atn_vs_cfp <- df_win |>
  filter(arm %in% c("atn", "cfp")) |>
  group_by(region, arm) |>
  summarise(cases = sum(clin_inc), .groups = "drop") |>
  pivot_wider(names_from = arm, values_from = cases) |>
  mutate(delta_clin = per1000_yr(atn - cfp))   # ATN minus Pyr-CFP, /1000/yr

scatter_df <- atn_vs_cfp |>
  left_join(eir_base, by = c(region = "name_1")) |>
  left_join(res_fut,  by = c(region = "name_1"))

# Symmetric fill limits so white sits exactly at 0 (no net difference).
fill_lim <- max(abs(scatter_df$delta_clin), na.rm = TRUE) * c(-1, 1)

p_scatter <- ggplot(scatter_df, aes(res_fut, eir, fill = delta_clin)) +
  geom_point(shape = 21, size = 6, colour = "grey25", stroke = 0.5) +
  ggrepel::geom_text_repel(aes(label = region), size = 3, colour = "grey20",
                           box.padding = 0.5, seed = 1, max.overlaps = Inf) +
  scale_y_log10() +
  # balance runs blue (low) -> white (0) -> red (high): negative delta
  # (ATN better) = blue, positive (ATN worse) = red. No direction flip.
  cmocean::scale_fill_cmocean(name = "balance", limits = fill_lim) +
  theme_minimal(base_size = 11) +
  labs(
    x = "Future pyrethroid resistance (mean over distribution years)",
    y = "Baseline EIR (log scale)",
    fill = "ATN - Pyr-CFP\nclinical cases\n(/1,000 / yr)",
    title = "Mali - ATN vs Pyr-CFP by resistance and transmission intensity",
    caption = paste0(
      "Each point = one admin-1 region. Fill = annual clinical incidence under ATN minus Pyr-CFP",
      " (per 1,000 / yr over the ", meta$n_future_years, "-yr window).\n",
      "Blue = ATN averts more than Pyr-CFP (better); red = ATN averts less (worse).",
      " x = arithmetic mean future resistance (", min(fut_years), "-", max(fut_years), ")."))

ggsave("dev/outputs/mali_resistance_eir_scatter.png", p_scatter,
       width = 9, height = 7, dpi = 150)
message("Saved: dev/outputs/mali_resistance_eir_scatter.png")

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

print(p_clin_series)
print(p_prev)
print(p_maps)
