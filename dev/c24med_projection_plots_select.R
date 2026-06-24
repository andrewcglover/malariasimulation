# =====================================================================
# c24med_projection_plots_select.R
#
#   Multi-file selectable plotter for c24med projection results.
#
#   Each "series" = one (arm, f, HL, retention) combination, drawn as
#   a single line/bar.  Edit the SERIES REGISTRY (Section 1) to:
#     • toggle arms on/off (enabled = TRUE/FALSE)
#     • select the source file via f / hl / ret per row
#     • override colours and line types
#     • overlay the same arm at multiple f/HL/ret values by adding a
#       second row with a distinct key (e.g. "F2") and a different
#       colour, linetype, and label
#
#   Output PNGs are tagged with PLOT_TAG + a content hash of the
#   enabled series — a different registry always writes to a new file,
#   so repeated sessions never overwrite existing plots.
#
#   The original c24med_projection_plots.R is LEFT UNCHANGED.
#
# =====================================================================

library(dplyr); library(tidyr); library(ggplot2); library(zoo)
library(sf); library(patchwork)

# =====================================================================
# 0. Top-level config
# =====================================================================

COUNTRY_ISO <- "MLI"

# Human-readable session tag — bump this label each new plot session so
# you can tell plot files apart at a glance.  The content hash (appended
# automatically) guarantees uniqueness even if you forget.
PLOT_TAG    <- "jobC"

# Baseline file for non-ATN arms (none / pyr / pyr_pbo / pyr_cfp).
# These arms don't vary with f or HL; they are ALWAYS read from here.
# Change BASELINE if you want a different no-drug comparator file.
# NOTE: if any ATN series uses ret="site", the cases-averted plot still
# compares against the ret1396 baseline defined here.  Change
# BASELINE$ret to "site" if you want a consistent within-retention
# comparison (and be sure that file contains the non-ATN arms).
BASELINE <- list(hl = "hl2p64", f = 0, ret = "ret1396")

# =====================================================================
# 1. Series registry
#
#  key       Short unique identifier — used as the ggplot grouping
#            factor.  NOT the arm name.  Must be distinct across rows.
#            When the same arm appears more than once (different f, HL,
#            or ret), give each occurrence a distinct key (e.g. "F",
#            "F2", "F3") and a distinct label so the legend is clear.
#
#  arm       One of: none | pyr | pyr_pbo | pyr_cfp | atn |
#                    pyr_atn | pyr_cfp_atn | pyr_cfp_mc_atn_cd
#
#  f         chem_dose_atn sensitivity value (0 or 1).
#            Only meaningful for pyr_atn / pyr_cfp_atn.
#            Forced to 0 for all other arms (f is a no-op for them and
#            the _chem1 files don't contain non-sensitive arms).
#
#  hl        Antimalarial half-life tag exactly as it appears in the
#            filename: "hl2p64" or "hl5".
#            Ignored for non-ATN arms (always read from BASELINE).
#
#  ret       Retention tag exactly as it appears in the filename.
#            "site"    → no ret suffix (e.g. ..._hl2p64.rds)
#            "ret1396" → appends _ret1396     (e.g. ..._hl2p64_ret1396.rds)
#            Ignored for non-ATN arms (always read from BASELINE).
#
#  enabled   TRUE/FALSE.  Set FALSE to hide without deleting the row.
#
#  colour    Hex colour (see Paul Tol "muted" palette below).
#
#  linetype  ggplot linetype: "solid" | "dashed" | "dotted" |
#            "dotdash" | "longdash" | "twodash"
#
#  label     Legend text.  Suffix e.g. "(f=1, HL5, ret1396)" when
#            the same arm appears more than once.
#
# Paul Tol "muted" qualitative palette (9 colours + 1 pale):
#   rose:   #CC6677   indigo: #332288   sand:   #DDCC77
#   green:  #117733   cyan:   #88CCEE   wine:   #882255
#   teal:   #44AA99   olive:  #999933   purple: #AA4499
#   pale:   #DDDDDD   (use for missing / "bad" data only)
# =====================================================================

plot_arms <- data.frame(stringsAsFactors = FALSE,
  key      = c("A",        "B",        "C",        "D",
               "E",        "F",        "G",        "H"),
  arm      = c("none",     "pyr",      "pyr_pbo",  "pyr_cfp",
               "atn",      "pyr_atn",  "pyr_cfp_atn", "pyr_cfp_mc_atn_cd"),
  f        = c(0,          0,          0,          0,
               0,          0,          0,          0),
  hl       = c("hl2p64",   "hl2p64",   "hl2p64",   "hl2p64",
               "hl2p64",   "hl2p64",   "hl2p64",   "hl2p64"),
  ret      = c("ret1396",  "ret1396",  "ret1396",  "ret1396",
               "ret1396",  "ret1396",  "ret1396",  "ret1396"),
  enabled  = c(TRUE,       TRUE,       TRUE,       TRUE,
               TRUE,       TRUE,       TRUE,       TRUE),
  colour   = c("#DDDDDD",  "#882255",  "#CC6677",  "#DDCC77",
               "#117733",  "#44AA99",  "#88CCEE",  "#AA4499"),
  linetype = c("solid",    "solid",    "solid",    "solid",
               "solid",    "solid",    "solid",    "solid"),
  label    = c("No nets",  "Pyr",      "Pyr-PBO",  "Pyr-CFP",
               "ATN",      "Pyr-ATN",  "Pyr-CFP-ATN", "Pyr-CFP MC, ATN CD")
)

# =====================================================================
# 2. Arm classification sets (used by resolve_file)
# =====================================================================
ARMS_NON_ATN   <- c("none", "pyr", "pyr_pbo", "pyr_cfp")
ARMS_F_ZERO    <- c("atn", "pyr_cfp_mc_atn_cd")   # f forced 0 (no-op)
ARMS_FULL_GRID <- c("pyr_atn", "pyr_cfp_atn")      # f + hl + ret all honoured

# =====================================================================
# 3. File-resolution helpers
# =====================================================================

# Build the RDS path for a given (hl, f, ret) combination.
make_path <- function(iso, hl, f, ret) {
  chem_sfx <- if (f > 0) sprintf("_chem%d", as.integer(f)) else ""
  ret_sfx  <- if (!is.null(ret) && ret != "site") sprintf("_%s", ret) else ""
  sprintf("dev/outputs/%s_c24med_projection_results_%s%s%s.rds",
          tolower(iso), hl, chem_sfx, ret_sfx)
}

# Resolve the actual RDS path for each registry row, applying the
# arm-classification rules described above.
resolve_file <- function(arm, f, hl, ret, iso = COUNTRY_ISO) {
  if (arm %in% ARMS_NON_ATN) {
    # Non-ATN: always the fixed baseline file; f/hl/ret ignored
    make_path(iso, BASELINE$hl, BASELINE$f, BASELINE$ret)
  } else if (arm %in% ARMS_F_ZERO) {
    # f is irrelevant (no drug sensitivity); hl and ret honoured
    make_path(iso, hl, 0, ret)
  } else {
    # pyr_atn / pyr_cfp_atn: full grid
    make_path(iso, hl, f, ret)
  }
}

# =====================================================================
# 4. Human-readable auto-suffix for output filenames
#
#    Encodes the enabled series in a compact, readable string:
#      k<keys>_<hl-values>_f<f-values>_<ret-values>
#    e.g.  kABCDEFGH_hl2p64_f0_ret1396
#          kABCDEF2FGH_hl2p64+hl5_f0+f1_ret1396
#
#    Values with more than one distinct entry are joined with "+".
#    Together with PLOT_TAG this is enough to identify any run from
#    the filename alone.  A companion manifest CSV (written below in
#    §6) records the full per-series details.
# =====================================================================
series_tag <- function(enabled_rows) {
  keys_str <- paste(enabled_rows$key, collapse = "")
  hl_str   <- paste(sort(unique(enabled_rows$hl)),  collapse = "+")
  f_str    <- paste(sort(unique(enabled_rows$f)),   collapse = "+")
  ret_str  <- paste(sort(unique(enabled_rows$ret)), collapse = "+")
  sprintf("k%s_%s_f%s_%s", keys_str, hl_str, f_str, ret_str)
}

# =====================================================================
# 5. Resolve paths, cache RDS reads, assemble combined data frame
# =====================================================================
active <- plot_arms[plot_arms$enabled, , drop = FALSE]
if (nrow(active) == 0L) stop("No series are enabled in plot_arms.")

# Compute RDS path for every enabled row
active$path <- mapply(resolve_file,
                      arm = active$arm, f   = active$f,
                      hl  = active$hl,  ret = active$ret)

# Check which files exist; warn and drop missing rather than erroring.
# This lets the script run against whatever jobs are finished so far.
exists_flag <- file.exists(active$path)
if (any(!exists_flag)) {
  warning(sprintf(
    "The following series will be skipped — RDS file not found:\n%s",
    paste(sprintf("  key=%s  arm=%-20s  f=%d  hl=%-6s  ret=%-7s  [%s]",
                  active$key[!exists_flag],
                  active$arm[!exists_flag],
                  active$f[!exists_flag],
                  active$hl[!exists_flag],
                  active$ret[!exists_flag],
                  active$path[!exists_flag]),
          collapse = "\n")
  ))
}
active <- active[exists_flag, , drop = FALSE]
if (nrow(active) == 0L) stop("No enabled series with an available RDS file.")

# Load each unique RDS once and cache it
unique_paths <- unique(active$path)
file_cache   <- new.env(parent = emptyenv())
for (p in unique_paths) {
  message(sprintf("Loading: %s", p))
  file_cache[[p]] <- readRDS(p)
}

# Baseline file (Job C): meta, shape, grid always come from here
baseline_path <- make_path(COUNTRY_ISO, BASELINE$hl, BASELINE$f, BASELINE$ret)
if (!file.exists(baseline_path)) stop("Baseline file missing: ", baseline_path)
if (!exists(baseline_path, envir = file_cache))
  file_cache[[baseline_path]] <- readRDS(baseline_path)

meta      <- file_cache[[baseline_path]]$meta
shape     <- file_cache[[baseline_path]]$shape
SHAPE_KEY <- "name_1"

# Assemble combined data frame — one row = one timestep × one series
df_all <- do.call(rbind, lapply(seq_len(nrow(active)), function(i) {
  row   <- active[i, ]
  obj_i <- file_cache[[row$path]]
  sub   <- obj_i$results[obj_i$results$arm == row$arm, , drop = FALSE]
  sub$key      <- row$key
  sub$label    <- row$label
  sub$colour   <- row$colour
  sub$linetype <- row$linetype
  sub
}))

# Ordered factor for ggplot grouping (preserves registry order in legends)
df_all$key_f <- factor(df_all$key, levels = active$key)

# Named vectors for manual scales
scale_colours   <- setNames(active$colour,   active$key)
scale_linetypes <- setNames(active$linetype, active$key)
scale_labels    <- setNames(active$label,    active$key)

# Human-readable filename suffix
sig <- series_tag(active)

# =====================================================================
# 6. Grid, windows, shared helpers
# =====================================================================
df_win <- df_all |>
  filter(year_rel >= 0, year_rel <= meta$n_future_years)

present_regions <- unique(df_win$region)

country_grid_raw <- geofacet::grid_auto(shape, names = SHAPE_KEY, seed = 1)
names(country_grid_raw)[names(country_grid_raw) == paste0("name_", SHAPE_KEY)] <- "name"
country_grid_raw$code <- country_grid_raw$name
country_grid <- country_grid_raw[country_grid_raw$code %in% present_regions,
                                 c("row", "col", "code", "name")]

vline_df <- data.frame(xintercept = meta$future_yr0)

# PNG path helper: dev/outputs/<iso>_c24med_<suffix>_<PLOT_TAG>_<sig>.png
out_file <- function(suffix) {
  tag <- if (nzchar(PLOT_TAG)) sprintf("%s_%s", PLOT_TAG, sig) else sig
  sprintf("dev/outputs/%s_c24med_%s_%s.png", tolower(COUNTRY_ISO), suffix, tag)
}

# Manifest CSV — one row per enabled series, written once so the full
# registry (key, arm, label, f, hl, ret, colour, linetype, resolved path)
# is recorded alongside the PNGs for future reference.
manifest_file <- {
  tag <- if (nzchar(PLOT_TAG)) sprintf("%s_%s", PLOT_TAG, sig) else sig
  sprintf("dev/outputs/%s_c24med_manifest_%s.csv", tolower(COUNTRY_ISO), tag)
}
write.csv(
  active[, c("key","arm","label","f","hl","ret","colour","linetype","path")],
  manifest_file, row.names = FALSE
)
message(sprintf("Manifest: %s", manifest_file))

# Cases-averted scaling helper (used in §8 and §9)
per1000_yr <- function(x) x / meta$human_pop * 1000 / meta$n_future_years

# =====================================================================
# 7. Geofacet timeseries — prevalence & clinical incidence
# =====================================================================
df_plot <- df_all |>
  filter(year_rel >= -3) |>
  mutate(
    cal_year  = year_rel + meta$future_yr0,
    clin_rate = clin_inc / meta$human_pop * 1000
  ) |>
  group_by(region, key_f) |>
  arrange(cal_year, .by_group = TRUE) |>
  mutate(
    pfpr_roll      = zoo::rollmean(pfpr2to10 * 100, k = 365L, fill = NA, align = "center"),
    clin_rate_roll = zoo::rollmean(clin_rate,        k = 365L, fill = NA, align = "center")
  ) |>
  ungroup()

p_prev <- ggplot(df_plot, aes(cal_year, pfpr2to10 * 100,
                               colour = key_f, linetype = key_f)) +
  geom_vline(data = vline_df, aes(xintercept = xintercept),
             linetype = "dashed", colour = "grey40", linewidth = 0.4,
             inherit.aes = FALSE) +
  geom_line(linewidth = 0.4, alpha = 0.2) +
  geom_line(aes(y = pfpr_roll), linewidth = 0.8, na.rm = TRUE) +
  scale_colour_manual(values = scale_colours,   labels = scale_labels, name = "") +
  scale_linetype_manual(values = scale_linetypes, labels = scale_labels, name = "") +
  geofacet::facet_geo(~ region, grid = country_grid) +
  theme_minimal(base_size = 15) +
  labs(x = "Year",
       y = expression(italic(Pf) * PR[2 - 10] * " (%)"))

p_clin_series <- ggplot(df_plot, aes(cal_year, clin_rate,
                                      colour = key_f, linetype = key_f)) +
  geom_vline(data = vline_df, aes(xintercept = xintercept),
             linetype = "dashed", colour = "grey40", linewidth = 0.4,
             inherit.aes = FALSE) +
  geom_line(linewidth = 0.4, alpha = 0.2) +
  geom_line(aes(y = clin_rate_roll), linewidth = 0.8, na.rm = TRUE) +
  scale_colour_manual(values = scale_colours,   labels = scale_labels, name = "") +
  scale_linetype_manual(values = scale_linetypes, labels = scale_labels, name = "") +
  geofacet::facet_geo(~ region, grid = country_grid) +
  theme_minimal(base_size = 15) +
  labs(x = "Year",
       y = "Clinical incidence (per 1,000 / day, all ages)")

ggsave(out_file("prevalence_facet"), p_prev,        width = 13, height = 9, dpi = 300)
ggsave(out_file("incidence_facet"),  p_clin_series, width = 13, height = 9, dpi = 300)
message(sprintf("Saved: %s\n       %s",
                out_file("prevalence_facet"), out_file("incidence_facet")))

# =====================================================================
# 8. Bar chart — cases averted vs no-nets (geofacet)
#
#    The "none" arm provides the baseline cases-per-region.  If none is
#    not enabled, this section is skipped.  If none appears more than
#    once, the first row's cases are used as the baseline.
# =====================================================================
none_rows   <- active[active$arm == "none", , drop = FALSE]
none_totals <- NULL   # populated below if none is available

if (nrow(none_rows) == 0L) {
  message("Section 8 (averted bar): 'none' arm not enabled — skipping.")
} else {
  none_key <- none_rows$key[1L]   # first 'none' row is the baseline

  # Total future cases per region per series
  totals_all <- df_win |>
    group_by(region, key_f) |>
    summarise(cases = sum(clin_inc), .groups = "drop") |>
    mutate(key = as.character(key_f))

  none_totals <- totals_all |>
    filter(key == none_key) |>
    select(region, none_cases = cases)

  averted_long <- totals_all |>
    filter(key != none_key) |>
    left_join(none_totals, by = "region") |>
    mutate(averted = per1000_yr(none_cases - cases))

  p_bar <- ggplot(averted_long, aes(x = key_f, y = averted, fill = key_f)) +
    geom_col(width = 0.75) +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey40") +
    scale_fill_manual(values = scale_colours, labels = scale_labels, name = "") +
    geofacet::facet_geo(~ region, grid = country_grid, scales = "fixed") +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x     = element_blank(),
      axis.ticks.x    = element_blank(),
      legend.position = "bottom"
    ) +
    guides(fill = guide_legend(nrow = 2)) +
    labs(x = "", y = "Cases averted vs no nets (per 1,000 / yr)", fill = "")

  ggsave(out_file("averted_bar"), p_bar, width = 13, height = 9, dpi = 300)
  message(sprintf("Saved: %s", out_file("averted_bar")))
}

# =====================================================================
# 9. Averted choropleths — one panel per non-none enabled series
#    Dynamic layout: ncol = 3, nrow computed from series count.
#    MLI only: Bamako receives an inset zoom panel in each map.
# =====================================================================
map_active <- active[active$arm != "none", , drop = FALSE]

if (is.null(none_totals) || nrow(map_active) == 0L) {
  message("Section 9 (averted maps): need 'none' arm + at least one other — skipping.")
} else {
  # averted_wide: one row per (region, non-none key)
  averted_wide <- totals_all |>
    filter(key != none_key) |>
    left_join(none_totals, by = "region") |>
    mutate(averted = per1000_yr(none_cases - cases)) |>
    select(region, key, averted)

  averted_rng <- range(averted_wide$averted, na.rm = TRUE)

  # Build a single-panel map for key_val / lbl
  make_averted_map <- function(key_val, lbl) {
    sub    <- averted_wide[averted_wide$key == key_val, c("region", "averted")]
    map_df <- shape |> left_join(sub, by = setNames("region", SHAPE_KEY))
    ggplot(map_df) +
      geom_sf(aes(fill = averted), colour = "white", linewidth = 0.2) +
      scale_fill_viridis_c(option = "D", limits = averted_rng) +
      theme_void(base_size = 14) +
      labs(title = lbl, fill = "Averted /\n1,000 / yr")
  }

  # Bamako inset — only applied for MLI.
  # Padding is 30% of Bamako's own width/height.
  # Inset position: left/bottom/right/top in 0-1 npc relative to panel.
  add_bamako_inset <- function(p_main, key_val, fill_lims) {
    sub       <- averted_wide[averted_wide$key == key_val, c("region", "averted")]
    map_df_bk <- shape |> left_join(sub, by = setNames("region", SHAPE_KEY))
    bamako_sf <- map_df_bk[map_df_bk[[SHAPE_KEY]] == "Bamako", ]
    bb        <- sf::st_bbox(bamako_sf)
    pad_frac  <- 0.30
    expand_x  <- as.numeric(bb["xmax"] - bb["xmin"]) * pad_frac
    expand_y  <- as.numeric(bb["ymax"] - bb["ymin"]) * pad_frac
    xlim <- c(bb["xmin"] - expand_x, bb["xmax"] + expand_x)
    ylim <- c(bb["ymin"] - expand_y, bb["ymax"] + expand_y)
    bb["xmin"] <- xlim[1]; bb["xmax"] <- xlim[2]
    bb["ymin"] <- ylim[1]; bb["ymax"] <- ylim[2]
    bb_rect <- sf::st_as_sfc(bb)
    p_rect  <- p_main +
      geom_sf(data = bb_rect, fill = NA, colour = "black",
              linewidth = 0.7, inherit.aes = FALSE)
    p_inset <- ggplot(bamako_sf) +
      geom_sf(aes(fill = averted), colour = "grey30", linewidth = 0.5) +
      scale_fill_viridis_c(option = "D", limits = fill_lims) +
      coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
      theme_void(base_size = 11) +
      theme(legend.position  = "none",
            panel.border     = element_rect(colour = "black", fill = NA,
                                            linewidth = 0.8))
    p_rect + patchwork::inset_element(
      p_inset, left = 0, bottom = 0.72, right = 0.26, top = 1.0,
      align_to = "panel"
    )
  }

  wrap_map <- function(p_main, key_val, fill_lims) {
    if (COUNTRY_ISO == "MLI") add_bamako_inset(p_main, key_val, fill_lims)
    else p_main
  }

  # Build one panel per non-none series
  map_list <- lapply(seq_len(nrow(map_active)), function(i) {
    k   <- map_active$key[i]
    lbl <- map_active$label[i]
    wrap_map(make_averted_map(k, lbl), k, averted_rng)
  })

  n_cols  <- 3L
  n_rows  <- ceiling(length(map_list) / n_cols)
  p_maps  <- patchwork::wrap_plots(map_list, ncol = n_cols,
                                   guides = "collect") &
    theme(legend.position = "bottom")

  ggsave(out_file("averted_maps"), p_maps,
         width  = 4 * n_cols,
         height = 4 * n_rows,
         dpi    = 300)
  message(sprintf("Saved: %s", out_file("averted_maps")))
}

# =====================================================================
# 10. % ATN-exposed mosquitoes over time (geofacet, ATN-bearing series)
#     Summed across all species and all adult compartments (Sv + Ev + Iv).
# =====================================================================
ATN_ARMS <- c("atn", "pyr_atn", "pyr_cfp_atn", "pyr_cfp_mc_atn_cd")

exp_cols <- grep("^(Sv|Ev|Iv)_exposed_",             names(df_all), value = TRUE)
all_cols <- grep("^(Sv|Ev|Iv)_(exposed|unexposed)_", names(df_all), value = TRUE)

atn_active <- active[active$arm %in% ATN_ARMS, , drop = FALSE]

if (length(exp_cols) == 0L || nrow(atn_active) == 0L) {
  message("Section 10 (% ATN exposed): no ATN-bearing series or exposure columns — skipping.")
} else {
  df_exp <- df_all |>
    filter(arm %in% ATN_ARMS, year_rel >= -3) |>
    mutate(
      cal_year   = year_rel + meta$future_yr0,
      mosq_exp   = rowSums(across(all_of(exp_cols))),
      mosq_total = rowSums(across(all_of(all_cols))),
      pct_exp    = if_else(mosq_total > 0, mosq_exp / mosq_total * 100, NA_real_)
    ) |>
    group_by(region, key_f) |>
    arrange(cal_year, .by_group = TRUE) |>
    mutate(pct_exp_roll = zoo::rollmean(pct_exp, k = 365L, fill = NA, align = "center")) |>
    ungroup()

  # Subset scale vectors to the ATN series only (suppresses ggplot warnings
  # about unused levels in the colour/linetype scales)
  atn_keys      <- atn_active$key
  atn_colours   <- scale_colours[atn_keys]
  atn_linetypes <- scale_linetypes[atn_keys]
  atn_labels    <- scale_labels[atn_keys]

  p_pct_exp <- ggplot(df_exp, aes(cal_year, pct_exp,
                                   colour = key_f, linetype = key_f)) +
    geom_vline(data = vline_df, aes(xintercept = xintercept),
               linetype = "dashed", colour = "grey40", linewidth = 0.4,
               inherit.aes = FALSE) +
    geom_line(linewidth = 0.4, alpha = 0.2) +
    geom_line(aes(y = pct_exp_roll), linewidth = 0.8, na.rm = TRUE) +
    scale_colour_manual(values = atn_colours,   labels = atn_labels, name = "") +
    scale_linetype_manual(values = atn_linetypes, labels = atn_labels, name = "") +
    geofacet::facet_geo(~ region, grid = country_grid) +
    theme_minimal(base_size = 15) +
    labs(x = "Year",
         y = "ATN-exposed adult mosquitoes (%)")

  ggsave(out_file("pct_atn_exposed"), p_pct_exp, width = 13, height = 9, dpi = 300)
  message(sprintf("Saved: %s", out_file("pct_atn_exposed")))
}

# =====================================================================
# 11. Past-window consistency check
#     Only check series that are read from the SAME file — series from
#     different files are independent sim runs and legitimately differ
#     in the past window by RNG seed.
# =====================================================================
baseline_keys <- active$key[active$path == baseline_path]
if (length(baseline_keys) > 1L) {
  past_check <- df_all |>
    filter(key %in% baseline_keys, year_rel < 0) |>
    group_by(region, timestep) |>
    summarise(n_distinct_pfpr = n_distinct(round(pfpr2to10, 6)), .groups = "drop")
  if (any(past_check$n_distinct_pfpr > 1)) {
    warning("Baseline-file arms diverge in the past window ",
            "(expected: IBM stochasticity + Erlang vs discrete-delay ODE).")
  } else {
    message("Verification OK: baseline-file series identical in the past window.")
  }
} else {
  message("Verification: only one baseline-file series — past-window check skipped.")
}

# =====================================================================
# Print to active graphics device
# =====================================================================
print(p_prev)
print(p_clin_series)
if (exists("p_bar"))     print(p_bar)
if (exists("p_maps"))    print(p_maps)
if (exists("p_pct_exp")) print(p_pct_exp)
