# =====================================================================
# mli_c24med_mosquito_density_plots.R
#   Geofacet (by admin-1 region) of total adult mosquito density expressed
#   as MOSQUITOES PER HUMAN (total M / human_pop), one line per arm.
#   Four separate plots: gambiae, arabiensis, funestus, and all-species total.
#   Log10 y-axis. Run from the fork root after the c24med projection results exist.
# =====================================================================

library(dplyr); library(tidyr); library(ggplot2); library(zoo); library(geofacet)

COUNTRY_ISO  <- "MLI"
RESULTS_FILE <- "dev/outputs/mli_c24med_projection_results.rds"

obj   <- readRDS(RESULTS_FILE)
df    <- obj$results
shape <- obj$shape
meta  <- obj$meta
hp    <- meta$human_pop
country_name <- if (!is.null(meta$country_name)) meta$country_name else COUNTRY_ISO
SHAPE_KEY <- "name_1"

# ---- arm labels / colours (7 arms in the c24med sweep) ----------------
arm_labels <- c(none        = "No future nets",
                pyr         = "Pyr-only",
                pyr_pbo     = "Pyr-PBO",
                pyr_cfp     = "Pyr-CFP",
                atn         = "ATN",
                pyr_atn     = "Pyr-ATN",
                pyr_cfp_atn = "Pyr-CFP-ATN")
arm_cols   <- c("No future nets" = "grey50",
                "Pyr-only"       = "#88CCEE",
                "Pyr-PBO"        = "#44AA99",
                "Pyr-CFP"        = "#009988",
                "ATN"            = "#EE7733",
                "Pyr-ATN"        = "#CC3311",
                "Pyr-CFP-ATN"    = "#882255")

# ---- total adult M per species, as ratio to humans --------------------
M_cols  <- function(sp) paste0(rep(c("Sv", "Ev", "Iv"), each = 2),
                               rep(c("_unexposed_", "_exposed_"), 3), sp, "_count")

dens <- df |>
  mutate(
    M_gambiae    = rowSums(across(all_of(M_cols("gambiae")))),
    M_arabiensis = rowSums(across(all_of(M_cols("arabiensis")))),
    M_funestus   = rowSums(across(all_of(M_cols("funestus")))),
    M_total      = M_gambiae + M_arabiensis + M_funestus
  ) |>
  mutate(across(starts_with("M_"), ~ . / hp)) |>          # ratio to humans
  mutate(arm_f    = factor(arm_labels[arm], levels = arm_labels),
         cal_year = year_rel + meta$future_yr0) |>
  filter(year_rel >= -3)

# ---- geofacet grid (regions present in the results) -------------------
grid_raw <- geofacet::grid_auto(shape, names = SHAPE_KEY, seed = 1)
names(grid_raw)[names(grid_raw) == paste0("name_", SHAPE_KEY)] <- "name"
grid_raw$code <- grid_raw$name
present  <- unique(dens$region)
country_grid <- grid_raw[grid_raw$code %in% present, c("row", "col", "code", "name")]

vline_df <- data.frame(xintercept = meta$future_yr0)

make_density_plot <- function(mcol, sp_label) {
  d <- dens |>
    group_by(region, arm_f) |>
    arrange(cal_year, .by_group = TRUE) |>
    mutate(M_roll = zoo::rollmean(.data[[mcol]], k = 365L, fill = NA, align = "center")) |>
    ungroup() |>
    filter(.data[[mcol]] > 0)   # log axis: drop non-positive (density is always > 0 anyway)

  ggplot(d, aes(cal_year, .data[[mcol]], colour = arm_f)) +
    geom_vline(data = vline_df, aes(xintercept = xintercept),
               linetype = "dashed", colour = "grey40", linewidth = 0.4) +
    geom_line(linewidth = 0.3, alpha = 0.20) +
    geom_line(aes(y = M_roll), linewidth = 0.7, na.rm = TRUE) +
    scale_colour_manual(values = arm_cols) +
    scale_y_log10() +
    annotation_logticks(sides = "l", colour = "grey70", linewidth = 0.25) +
    geofacet::facet_geo(~ region, grid = country_grid) +
    theme_minimal(base_size = 11) +
    labs(x = "Year", y = "Adult mosquitoes per human (log scale)", colour = "",
         title = paste0(country_name, " — ", sp_label,
                        " adult density (per human) by region and future net scenario"),
         caption = paste0("Dashed line = future distribution start (", meta$future_yr0,
                          "). Thick line = 365-day rolling mean. Log10 y-axis."))
}

plots <- list(
  gambiae    = make_density_plot("M_gambiae",    "An. gambiae"),
  arabiensis = make_density_plot("M_arabiensis", "An. arabiensis"),
  funestus   = make_density_plot("M_funestus",   "An. funestus"),
  total      = make_density_plot("M_total",      "All-species total")
)

for (nm in names(plots)) {
  out <- sprintf("dev/outputs/mli_c24med_mosquito_density_%s.png", nm)
  ggsave(out, plots[[nm]], width = 11, height = 8, dpi = 150)
  message("Saved: ", out)
}
