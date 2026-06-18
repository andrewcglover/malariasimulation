# dev/mali_single_region_test.R
#
# Smoke test: run all 4 arms sequentially for the highest-EIR Mali admin-1
# region, then print a summary table.
#
# Usage (from the fork root, after devtools::load_all()):
#   source("dev/mali_single_region_test.R")

# Source the run script for its functions + data; skip the cluster launch.
options(mali_test_mode = TRUE)
source("dev/mali_projection_run.R")
options(mali_test_mode = NULL)

future_start_year <- 2025

# ── Pick highest-EIR region ───────────────────────────────────────────────────
eir_by_region <- vapply(regions, function(rg) {
  site_row <- site_obj$sites[site_obj$sites$name_1 == rg, , drop = FALSE]
  ms       <- site::subset_site(site_obj, site_row)
  ms$eir$eir
}, numeric(1))

test_region <- names(which.max(eir_by_region))
message(sprintf("\nHighest-EIR region: %s  (EIR = %.1f)\n", test_region, max(eir_by_region)))

# ── Run all 4 arms for that region ───────────────────────────────────────────
results <- lapply(arms, function(arm) {
  message(sprintf("  arm = %-8s ...", arm), appendLF = FALSE)
  t0  <- proc.time()[["elapsed"]]
  out <- run_one(list(region = test_region, arm = arm))
  dt  <- proc.time()[["elapsed"]] - t0
  message(sprintf(" done (%.0f s)", dt))
  out
})

df <- dplyr::bind_rows(results)

max_rel_yr <- max(df$year_rel)

# ── Summary: mean future PfPR and clin incidence per arm ─────────────────────
future_df <- df[df$year_rel >= 0, ]
summary_tbl <- future_df |>
  dplyr::group_by(arm) |>
  dplyr::summarise(
    mean_pfpr2to10  = round(mean(pfpr2to10,  na.rm = TRUE), 4),
    mean_clin_inc   = round(sum(clin_inc,   na.rm = TRUE) / max_rel_yr, 4),
    .groups = "drop"
  )

message("\nFuture-period summary (years 0–6 relative to ", future_yr0, "):")
print(as.data.frame(summary_tbl))

# ── Plots ─────────────────────────────────────────────────────────────────────
library(ggplot2)

arm_labels <- c(none    = "No future nets",
                cfp     = "Future Pyr-CFP",
                atn     = "Future ATN",
                pyr_atn = "Future Pyr-ATN")
arm_cols   <- c("No future nets" = "grey50",
                "Future Pyr-CFP" = "#009988",
                "Future ATN"     = "#EE7733",
                "Future Pyr-ATN" = "#CC3311")

df_plot <- df |>
  dplyr::filter(year_rel >= -3) |>   # show 5 years of history + future window
  dplyr::mutate(arm_f = factor(arm_labels[arm], levels = arm_labels))

vline_df <- data.frame(xintercept = future_start_year)

# 0. PfPR 2-10 over time
p_prev <- ggplot(df_plot, aes(year_rel + future_start_year, pfpr2to10 * 100, colour = arm_f)) +
  geom_vline(data = vline_df, aes(xintercept = xintercept),
             linetype = "dashed", colour = "grey40", linewidth = 0.4) +
  geom_line(linewidth = 0.7) +
  scale_colour_manual(values = arm_cols) +
  theme_minimal(base_size = 12) +
  labs(x = "Years relative to first future distribution",
       y = expression(italic(Pf) * PR[2 - 10] * " (%)"),
       colour = "",
       title  = sprintf("%s", test_region))

print(p_prev)

# 1. Net use over time
p_use <- ggplot(df_plot, aes(year_rel + future_start_year, n_use_net / human_pop * 100, colour = arm_f)) +
  geom_vline(data = vline_df, aes(xintercept = xintercept),
             linetype = "dashed", colour = "grey40", linewidth = 0.4) +
  geom_line(linewidth = 0.7) +
  scale_colour_manual(values = arm_cols) +
  theme_minimal(base_size = 12) +
  labs(x = "Years relative to first future distribution",
       y = "Net use (%)",
       colour = "",
       title  = sprintf("%s", test_region))

print(p_use)

# 2. Clinical incidence over time (all ages, per 1000 pop / yr)
df_plot2 <- df_plot |>
  dplyr::mutate(clin_inc_rate = clin_inc / human_pop * 1000)

p_inc <- ggplot(df_plot2, aes(year_rel + future_start_year, clin_inc_rate, colour = arm_f)) +
  geom_vline(data = vline_df, aes(xintercept = xintercept),
             linetype = "dashed", colour = "grey40", linewidth = 0.4) +
  geom_line(linewidth = 0.7) +
  scale_colour_manual(values = arm_cols) +
  theme_minimal(base_size = 12) +
  labs(x = "Years relative to first future distribution",
       y = "Clinical incidence (per 1,000 / yr, all ages)",
       colour = "",
       title  = sprintf("%s", test_region))

print(p_inc)

# 3. Cases averted vs 'none' arm over the future window (bar chart)
none_cases <- summary_tbl$mean_clin_inc[summary_tbl$arm == "none"]
averted_tbl <- summary_tbl |>
  dplyr::filter(arm != "none") |>
  dplyr::mutate(
    arm_f      = factor(arm_labels[arm], levels = arm_labels),
    averted_rate = (none_cases - mean_clin_inc) / human_pop * 1000
  )

p_avert <- ggplot(averted_tbl, aes(arm_f, averted_rate, fill = arm_f)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  scale_fill_manual(values = arm_cols) +
  theme_minimal(base_size = 12) +
  labs(x = "",
       y = "Clinical cases averted vs no-nets (per 1,000 / yr)",
       title  = sprintf("%s", test_region))

print(p_avert)

# 4. EIR
p_EIR_arabiensis <- ggplot(df_plot, aes(year_rel + future_start_year, EIR_arabiensis, colour = arm_f)) +
  geom_vline(data = vline_df, aes(xintercept = xintercept),
             linetype = "dashed", colour = "grey40", linewidth = 0.4) +
  geom_line(linewidth = 0.7) +
  scale_colour_manual(values = arm_cols) +
  theme_minimal(base_size = 12) +
  labs(x = "Years relative to first future distribution",
       y = "EIR arabiensis",
       colour = "",
       title  = sprintf("%s", test_region))

print(p_EIR_arabiensis)

p_EIR_funestus <- ggplot(df_plot, aes(year_rel + future_start_year, EIR_funestus, colour = arm_f)) +
  geom_vline(data = vline_df, aes(xintercept = xintercept),
             linetype = "dashed", colour = "grey40", linewidth = 0.4) +
  geom_line(linewidth = 0.7) +
  scale_colour_manual(values = arm_cols) +
  theme_minimal(base_size = 12) +
  labs(x = "Years relative to first future distribution",
       y = "EIR funestus",
       colour = "",
       title  = sprintf("%s", test_region))

print(p_EIR_funestus)

p_EIR_gambiae <- ggplot(df_plot, aes(year_rel + future_start_year, EIR_gambiae, colour = arm_f)) +
  geom_vline(data = vline_df, aes(xintercept = xintercept),
             linetype = "dashed", colour = "grey40", linewidth = 0.4) +
  geom_line(linewidth = 0.7) +
  scale_colour_manual(values = arm_cols) +
  theme_minimal(base_size = 12) +
  labs(x = "Year",
       y = "EIR gambiae",
       colour = "",
       title  = sprintf("%s", test_region))

print(p_EIR_gambiae)

# 4b. EIR by species, combined plot

eir_plot <- df_plot |>
  dplyr::select(
    year_rel, arm_f,
    dplyr::matches("^EIR_(gambiae|arabiensis|funestus)$")
  ) |>
  tidyr::pivot_longer(
    cols = dplyr::matches("^EIR_(gambiae|arabiensis|funestus)$"),
    names_to = "species",
    names_prefix = "EIR_",
    values_to = "EIR"
  ) |>
  dplyr::mutate(
    species = factor(
      species,
      levels = c("gambiae", "arabiensis", "funestus"),
      labels = c("gambiae", "arabiensis", "funestus")
    )
  )

p_eir_all <- ggplot(
  eir_plot,
  aes(
    x = year_rel + future_start_year,
    y = EIR,
    colour = arm_f,
    group = arm_f
  )
) +
  geom_vline(
    data = vline_df,
    aes(xintercept = xintercept),
    linetype = "dashed",
    colour = "grey40",
    linewidth = 0.4
  ) +
  geom_line(linewidth = 0.7) +
  facet_wrap(~ species, nrow = 1) +
  scale_colour_manual(values = arm_cols) +
  theme_minimal(base_size = 12) +
  labs(
    x = "Years relative to first future distribution",
    y = "EIR",
    colour = "",
    title = sprintf("%s", test_region)
  )

print(p_eir_all)

# 5. Vector counts by species and compartment
#    Solid = all mosquitoes (ATN-exposed + unexposed)
#    Dashed = ATN-exposed only
#    Difference between solid and dashed = ATN-unexposed implicitly

species_vec <- c("gambiae", "arabiensis", "funestus")
species_labs <- c(
  gambiae    = "gambiae",
  arabiensis = "arabiensis",
  funestus   = "funestus"
)

df_vec <- df |>
  dplyr::filter(year_rel >= -3) |>
  dplyr::mutate(
    row_id = dplyr::row_number(),
    human_pop = n_age_0_1824 + n_age_1825_5474 + n_age_5475_36499,
    arm_f = factor(arm_labels[arm], levels = arm_labels)
  )

vector_long <- df_vec |>
  dplyr::select(
    row_id, year_rel, arm, arm_f, region, human_pop,
    dplyr::matches("^(Sv|Ev|Iv)_(unexposed|exposed)_(gambiae|arabiensis|funestus)_count$")
  ) |>
  tidyr::pivot_longer(
    cols = dplyr::matches("^(Sv|Ev|Iv)_(unexposed|exposed)_(gambiae|arabiensis|funestus)_count$"),
    names_to = c("compartment", "exposure", "species"),
    names_pattern = "^(Sv|Ev|Iv)_(unexposed|exposed)_(gambiae|arabiensis|funestus)_count$",
    values_to = "count"
  )

compartment_dat <- vector_long |>
  dplyr::group_by(row_id, year_rel, arm_f, species, compartment) |>
  dplyr::summarise(
    all_count = sum(count, na.rm = TRUE),
    atn_exposed_count = sum(count[exposure == "exposed"], na.rm = TRUE),
    human_pop = dplyr::first(human_pop),
    .groups = "drop"
  )

total_dat <- vector_long |>
  dplyr::group_by(row_id, year_rel, arm_f, species) |>
  dplyr::summarise(
    compartment = "total",
    all_count = sum(count, na.rm = TRUE),
    atn_exposed_count = sum(count[exposure == "exposed"], na.rm = TRUE),
    human_pop = dplyr::first(human_pop),
    .groups = "drop"
  )

plot_dat <- dplyr::bind_rows(compartment_dat, total_dat) |>
  tidyr::pivot_longer(
    cols = c(all_count, atn_exposed_count),
    names_to = "line",
    values_to = "count"
  ) |>
  dplyr::mutate(
    line = dplyr::recode(
      line,
      all_count = "All mosquitoes",
      atn_exposed_count = "ATN-exposed only"
    ),
    compartment = factor(compartment, levels = c("Sv", "Ev", "Iv", "total")),
    species = factor(species, levels = c("gambiae", "arabiensis", "funestus"))
  ) |>
  dplyr::group_by(year_rel, arm_f, species, compartment, line) |>
  dplyr::summarise(
    count = sum(count, na.rm = TRUE),
    human_pop = sum(human_pop, na.rm = TRUE),
    mosquitoes_per_human = count / human_pop,
    .groups = "drop"
  )

# Species-specific plots (4 facets each: Sv, Ev, Iv, total)
plot_one_species <- function(sp) {
  ggplot(
    dplyr::filter(plot_dat, species == sp),
    aes(
      x = year_rel + future_start_year,
      y = mosquitoes_per_human,
      colour = arm_f,
      linetype = line,
      group = interaction(arm_f, line)
    )
  ) +
    geom_vline(
      data = vline_df,
      aes(xintercept = xintercept),
      linetype = "dashed",
      colour = "grey40",
      linewidth = 0.4
    ) +
    geom_line(linewidth = 0.7) +
    facet_wrap(~ compartment, scales = "free_y", ncol = 2) +
    scale_colour_manual(values = arm_cols) +
    scale_linetype_manual(
      values = c(
        "All mosquitoes" = "solid",
        "ATN-exposed only" = "dashed"
      )
    ) +
    theme_minimal(base_size = 12) +
    labs(
      x = "Years relative to first future distribution",
      y = "Mosquitoes per human",
      colour = "",
      linetype = "",
      title = sprintf("%s: %s", test_region, species_labs[[sp]])
    )
}

p_vec_gambiae    <- plot_one_species("gambiae")
p_vec_arabiensis <- plot_one_species("arabiensis")
p_vec_funestus   <- plot_one_species("funestus")

print(p_vec_gambiae)
print(p_vec_arabiensis)
print(p_vec_funestus)

# Combined facet plot: rows = species, cols = compartment
p_vec_all <- ggplot(
  plot_dat,
  aes(
    x = year_rel + future_start_year,
    y = mosquitoes_per_human,
    colour = arm_f,
    linetype = line,
    group = interaction(arm_f, line)
  )
) +
  geom_vline(
    data = vline_df,
    aes(xintercept = xintercept),
    linetype = "dashed",
    colour = "grey40",
    linewidth = 0.4
  ) +
  geom_line(linewidth = 0.7) +
  facet_grid(compartment ~ species, scales = "free_y") +
  scale_colour_manual(values = arm_cols) +
  scale_linetype_manual(
    values = c(
      "All mosquitoes" = "solid",
      "ATN-exposed only" = "dashed"
    )
  ) +
  theme_minimal(base_size = 12) +
  labs(
    x = "Year",
    y = "Adult female mosquitoes per human",
    colour = "",
    linetype = "",
    title = sprintf("%s", test_region)
  )

print(p_vec_all)