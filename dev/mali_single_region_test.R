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
library(zoo)

arm_labels <- c(none    = "No future nets",
                cfp     = "Future Pyr-CFP",
                atn     = "Future ATN",
                pyr_atn = "Future Pyr-ATN")
arm_cols   <- c("No future nets" = "grey50",
                "Future Pyr-CFP" = "#009988",
                "Future ATN"     = "#EE7733",
                "Future Pyr-ATN" = "#CC3311")

df_plot <- df |>
  dplyr::filter(year_rel >= -3) |>   # show 3 years of history + future window
  dplyr::mutate(arm_f = factor(arm_labels[arm], levels = arm_labels))

vline_df <- data.frame(xintercept = future_start_year)

# ── Plot output directory ─────────────────────────────────────────────────────
plot_dir <- file.path("dev", "outputs", "single_region_plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

save_plot <- function(p, name, width = 8, height = 5) {
  ggsave(
    file.path(plot_dir, sprintf("%s_%s.png", test_region, name)),
    p, width = width, height = height, dpi = 150
  )
}

# ── Rolling 365-day means (centre-aligned) for df_plot ───────────────────────
df_plot <- df_plot |>
  dplyr::group_by(arm_f) |>
  dplyr::arrange(year_rel, .by_group = TRUE) |>
  dplyr::mutate(
    pfpr_pct_roll          = zoo::rollmean(pfpr2to10 * 100,              k = 365L, fill = NA, align = "center"),
    net_use_pct_roll       = zoo::rollmean(n_use_net / human_pop * 100,  k = 365L, fill = NA, align = "center"),
    clin_rate_roll         = zoo::rollmean(clin_inc  / human_pop * 1000, k = 365L, fill = NA, align = "center"),
    EIR_gambiae_pp_roll    = zoo::rollmean(EIR_gambiae_pp,               k = 365L, fill = NA, align = "center"),
    EIR_arabiensis_pp_roll = zoo::rollmean(EIR_arabiensis_pp,            k = 365L, fill = NA, align = "center"),
    EIR_funestus_pp_roll   = zoo::rollmean(EIR_funestus_pp,              k = 365L, fill = NA, align = "center")
  ) |>
  dplyr::ungroup()

# 0. PfPR 2-10 over time
p_prev <- ggplot(df_plot, aes(year_rel + future_start_year, pfpr2to10 * 100, colour = arm_f)) +
  geom_vline(data = vline_df, aes(xintercept = xintercept),
             linetype = "dashed", colour = "grey40", linewidth = 0.4) +
  geom_line(linewidth = 0.4, alpha = 0.2) +
  geom_line(aes(y = pfpr_pct_roll), linewidth = 0.8, na.rm = TRUE) +
  scale_colour_manual(values = arm_cols) +
  theme_minimal(base_size = 12) +
  labs(x = "Year",
       y = expression(italic(Pf) * PR[2 - 10] * " (%)"),
       colour = "",
       title  = sprintf("%s", test_region))

print(p_prev)
save_plot(p_prev, "prev")

# 1. Net use over time
p_use <- ggplot(df_plot, aes(year_rel + future_start_year, n_use_net / human_pop * 100, colour = arm_f)) +
  geom_vline(data = vline_df, aes(xintercept = xintercept),
             linetype = "dashed", colour = "grey40", linewidth = 0.4) +
  geom_line(linewidth = 0.4, alpha = 0.2) +
  geom_line(aes(y = net_use_pct_roll), linewidth = 0.8, na.rm = TRUE) +
  scale_colour_manual(values = arm_cols) +
  theme_minimal(base_size = 12) +
  labs(x = "Year",
       y = "Net use (%)",
       colour = "",
       title  = sprintf("%s", test_region))

print(p_use)
save_plot(p_use, "net_use")

# 2. Clinical incidence over time (all ages, per 1000 pop / yr)
df_plot2 <- df_plot |>
  dplyr::mutate(clin_inc_rate = clin_inc / human_pop * 1000)
# clin_rate_roll is inherited from df_plot

p_inc <- ggplot(df_plot2, aes(year_rel + future_start_year, clin_inc_rate, colour = arm_f)) +
  geom_vline(data = vline_df, aes(xintercept = xintercept),
             linetype = "dashed", colour = "grey40", linewidth = 0.4) +
  geom_line(linewidth = 0.4, alpha = 0.2) +
  geom_line(aes(y = clin_rate_roll), linewidth = 0.8, na.rm = TRUE) +
  scale_colour_manual(values = arm_cols) +
  theme_minimal(base_size = 12) +
  labs(x = "Year",
       y = "Clinical incidence (per 1,000 / yr, all ages)",
       colour = "",
       title  = sprintf("%s", test_region))

print(p_inc)
save_plot(p_inc, "clin_inc")

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
save_plot(p_avert, "cases_averted")

# 4. EIR (bites / person / day)
p_EIR_arabiensis <- ggplot(df_plot, aes(year_rel + future_start_year, EIR_arabiensis_pp, colour = arm_f)) +
  geom_vline(data = vline_df, aes(xintercept = xintercept),
             linetype = "dashed", colour = "grey40", linewidth = 0.4) +
  geom_line(linewidth = 0.4, alpha = 0.2) +
  geom_line(aes(y = EIR_arabiensis_pp_roll), linewidth = 0.8, na.rm = TRUE) +
  scale_colour_manual(values = arm_cols) +
  theme_minimal(base_size = 12) +
  labs(x = "Year",
       y = "EIR arabiensis (bites/person/day)",
       colour = "",
       title  = sprintf("%s", test_region))

print(p_EIR_arabiensis)
save_plot(p_EIR_arabiensis, "eir_arabiensis")

p_EIR_funestus <- ggplot(df_plot, aes(year_rel + future_start_year, EIR_funestus_pp, colour = arm_f)) +
  geom_vline(data = vline_df, aes(xintercept = xintercept),
             linetype = "dashed", colour = "grey40", linewidth = 0.4) +
  geom_line(linewidth = 0.4, alpha = 0.2) +
  geom_line(aes(y = EIR_funestus_pp_roll), linewidth = 0.8, na.rm = TRUE) +
  scale_colour_manual(values = arm_cols) +
  theme_minimal(base_size = 12) +
  labs(x = "Year",
       y = "EIR funestus (bites/person/day)",
       colour = "",
       title  = sprintf("%s", test_region))

print(p_EIR_funestus)
save_plot(p_EIR_funestus, "eir_funestus")

p_EIR_gambiae <- ggplot(df_plot, aes(year_rel + future_start_year, EIR_gambiae_pp, colour = arm_f)) +
  geom_vline(data = vline_df, aes(xintercept = xintercept),
             linetype = "dashed", colour = "grey40", linewidth = 0.4) +
  geom_line(linewidth = 0.4, alpha = 0.2) +
  geom_line(aes(y = EIR_gambiae_pp_roll), linewidth = 0.8, na.rm = TRUE) +
  scale_colour_manual(values = arm_cols) +
  theme_minimal(base_size = 12) +
  labs(x = "Year",
       y = "EIR gambiae (bites/person/day)",
       colour = "",
       title  = sprintf("%s", test_region))

print(p_EIR_gambiae)
save_plot(p_EIR_gambiae, "eir_gambiae")

# 4b. EIR by species, combined plot (bites / person / day)
eir_plot <- df_plot |>
  dplyr::select(
    year_rel, arm_f,
    dplyr::matches("^EIR_(gambiae|arabiensis|funestus)_pp$")
  ) |>
  tidyr::pivot_longer(
    cols = dplyr::matches("^EIR_(gambiae|arabiensis|funestus)_pp$"),
    names_to = "species",
    names_pattern = "^EIR_(gambiae|arabiensis|funestus)_pp$",
    values_to = "EIR_pp"
  ) |>
  dplyr::mutate(
    species = factor(
      species,
      levels = c("gambiae", "arabiensis", "funestus"),
      labels = c("gambiae", "arabiensis", "funestus")
    )
  ) |>
  dplyr::group_by(arm_f, species) |>
  dplyr::arrange(year_rel, .by_group = TRUE) |>
  dplyr::mutate(EIR_pp_roll = zoo::rollmean(EIR_pp, k = 365L, fill = NA, align = "center")) |>
  dplyr::ungroup()

p_eir_all <- ggplot(
  eir_plot,
  aes(
    x = year_rel + future_start_year,
    y = EIR_pp,
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
  geom_line(linewidth = 0.4, alpha = 0.2) +
  geom_line(aes(y = EIR_pp_roll), linewidth = 0.8, na.rm = TRUE) +
  facet_wrap(~ species, nrow = 1) +
  scale_colour_manual(values = arm_cols) +
  theme_minimal(base_size = 12) +
  labs(
    x = "Year",
    y = "EIR (bites/person/day)",
    colour = "",
    title = sprintf("%s", test_region)
  )

print(p_eir_all)
save_plot(p_eir_all, "eir_all_species", width = 12, height = 5)

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
  ) |>
  dplyr::group_by(arm_f, species, compartment, line) |>
  dplyr::arrange(year_rel, .by_group = TRUE) |>
  dplyr::mutate(mosq_roll = zoo::rollmean(mosquitoes_per_human, k = 365L, fill = NA, align = "center")) |>
  dplyr::ungroup()

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
    geom_line(linewidth = 0.4, alpha = 0.2) +
    geom_line(aes(y = mosq_roll), linewidth = 0.8, na.rm = TRUE) +
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
      x = "Year",
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

save_plot(p_vec_gambiae,    "vec_gambiae",    width = 8, height = 7)
save_plot(p_vec_arabiensis, "vec_arabiensis", width = 8, height = 7)
save_plot(p_vec_funestus,   "vec_funestus",   width = 8, height = 7)

# Combined facet plot: rows = compartment, cols = species
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
  geom_line(linewidth = 0.4, alpha = 0.2) +
  geom_line(aes(y = mosq_roll), linewidth = 0.8, na.rm = TRUE) +
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
save_plot(p_vec_all, "vec_all", width = 14, height = 10)

# 5b. Total mosquitoes only, by species (compartment ~ species grid, total row only)
p_vec_total <- ggplot(
  dplyr::filter(plot_dat, compartment == "total"),
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
  geom_line(linewidth = 0.4, alpha = 0.2) +
  geom_line(aes(y = mosq_roll), linewidth = 0.8, na.rm = TRUE) +
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

print(p_vec_total)
save_plot(p_vec_total, "vec_total", width = 12, height = 4)

# 6. Percentage ATN-exposed by species and compartment

plot_dat_exp <- dplyr::bind_rows(compartment_dat, total_dat) |>
  dplyr::mutate(
    pct_atn_exposed = dplyr::if_else(
      all_count > 0,
      atn_exposed_count / all_count * 100,
      NA_real_
    ),
    compartment = factor(compartment, levels = c("Sv", "Ev", "Iv", "total")),
    species = factor(species, levels = c("gambiae", "arabiensis", "funestus"))
  ) |>
  dplyr::group_by(year_rel, arm_f, species, compartment) |>
  dplyr::summarise(
    atn_exposed_count = sum(atn_exposed_count, na.rm = TRUE),
    all_count = sum(all_count, na.rm = TRUE),
    pct_atn_exposed = dplyr::if_else(
      all_count > 0,
      atn_exposed_count / all_count * 100,
      NA_real_
    ),
    .groups = "drop"
  ) |>
  dplyr::group_by(arm_f, species, compartment) |>
  dplyr::arrange(year_rel, .by_group = TRUE) |>
  dplyr::mutate(
    pct_roll = zoo::rollmean(pct_atn_exposed, k = 365L, fill = NA, align = "center"),
    pct_roll = dplyr::if_else(year_rel >= 0.5, pct_roll, NA_real_)
  ) |>
  dplyr::ungroup()

plot_pct_exposed_one_species <- function(sp) {
  ggplot(
    dplyr::filter(plot_dat_exp, species == sp,
                  arm_f %in% c("Future ATN", "Future Pyr-ATN")),
    aes(
      x = year_rel + future_start_year,
      y = pct_atn_exposed,
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
    geom_line(linewidth = 0.4, alpha = 0.2) +
    geom_line(aes(y = pct_roll), linewidth = 0.8, na.rm = TRUE) +
    facet_wrap(~ compartment, scales = "free_y", ncol = 2) +
    scale_colour_manual(values = arm_cols) +
    theme_minimal(base_size = 12) +
    labs(
      x = "Year",
      y = "ATN-exposed mosquitoes (%)",
      colour = "",
      title = sprintf("%s: %s", test_region, species_labs[[sp]])
    )
}

p_vec_exp_gambiae    <- plot_pct_exposed_one_species("gambiae")
p_vec_exp_arabiensis <- plot_pct_exposed_one_species("arabiensis")
p_vec_exp_funestus   <- plot_pct_exposed_one_species("funestus")

print(p_vec_exp_gambiae)
print(p_vec_exp_arabiensis)
print(p_vec_exp_funestus)

save_plot(p_vec_exp_gambiae,    "vec_exp_gambiae",    width = 8, height = 7)
save_plot(p_vec_exp_arabiensis, "vec_exp_arabiensis", width = 8, height = 7)
save_plot(p_vec_exp_funestus,   "vec_exp_funestus",   width = 8, height = 7)

p_vec_exp_all <- ggplot(
  dplyr::filter(plot_dat_exp, arm_f %in% c("Future ATN", "Future Pyr-ATN")),
  aes(
    x = year_rel + future_start_year,
    y = pct_atn_exposed,
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
  geom_line(linewidth = 0.4, alpha = 0.2) +
  geom_line(aes(y = pct_roll), linewidth = 0.8, na.rm = TRUE) +
  facet_grid(compartment ~ species, scales = "free_y") +
  scale_colour_manual(values = arm_cols) +
  theme_minimal(base_size = 12) +
  labs(
    x = "Year",
    y = "ATN-exposed mosquitoes (%)",
    colour = "",
    title = sprintf("%s", test_region)
  )

print(p_vec_exp_all)
save_plot(p_vec_exp_all, "vec_exp_all", width = 14, height = 10)
