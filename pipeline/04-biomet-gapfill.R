# Gap-fill biomet variables ----------------------------------------------------

# Purpose: 

# Load the required packages
library("REddyProc")
library("openeddy")
library("yaml")
library("glue")
library("tidyverse")

# Load custom functions
source("R/pipeline-funs.R")

## Initialize script settings & documentation ----------------------------------

if (interactive()) {
  site <- "JLA" # three-letter site code
  years <- 2019:2025 # years to process
}
# years <- 2019:2025

output_dirs <- file.path("output", site, years) # where to write output
report_dir <- file.path("reports", site) # where to save reports

# Load metadata file
md <- read_yaml(file.path("data", site, "metadata.yml"))

# Load variable info file
vars <- read_yaml("config/var-config.yml")

# Load configuration file
cfg <- read_yaml("config/config.yml")

# Set configuration options
na_frac <- 0.5 # min % available points for valid de-biasing coefs
max_interp <- 3 # max gap length for linear interpolation
# max_kalman <- 48 # max gap length for kalman filter interpolation
precip_thr <- 0.1 # >0 since ERA interpolation can create tiny values

# Biomet vars to be filled
gf_vars <- cfg$met_vars
soil_vars <- cfg$biomet_gapfill$soil_vars

## Read required data files ----------------------------------------------------

# Quality-controlled data
data_files <- glue("eddypro_{site}-{years}_fluxnet_adv_qaqc.csv")

# data <- read_csv(file.path(output_dirs, data_files))
data <- output_dirs |>
  file.path(data_files) |>
  map(read_csv) |>
  bind_rows()

# Alternative site data (quality-controlled)
sites <- cfg$sites
alt_sites <- sites[sites != site]
alt_dirs <- map(alt_sites, \(x) file.path("output", x, years))
alt_files <- map(
  alt_sites, \(x) glue("eddypro_{x}-{years}_fluxnet_adv_qaqc.csv")
)
# alt_data <- alt_dirs |> 
#   map2(alt_files, \(x, y) file.path(x, y)) |>
#   map(read_csv) |>
#   set_names(alt_sites)
alt_data <- alt_dirs |> 
  map2(alt_files, \(x, y) file.path(x, y)) |>
  map(as.list) |>
  map_depth(2, read_csv) |>
  map(bind_rows) |>
  set_names(alt_sites)

# ERA data
era <- read_csv(glue("output/ERA5/era5-sl-{years}-hh.csv"))

## Clean biomet vars with QC flags ---------------------------------------------

cleaned <- gf_vars |>
  map2(
    glue("QC_{gf_vars}"),
    \(x, y) transmute(data, "{x}" := if_else(.data[[y]] == 2, NA, .data[[x]]))
  ) |>
  bind_cols() |>
  mutate(across(where(is.logical), as.numeric)) |>
  mutate(TIMESTAMP = data$TIMESTAMP, .before = 1)

# Data coverage plot
# TODO: compress these somehow; too much vector going on
coverage_plot <- cleaned |>
  pivot_longer(-TIMESTAMP) |>
  # drop_na() |>
  mutate(coverage = as.integer(!is.na(value))) |>
  nest_by(year = year(TIMESTAMP)) |>
  mutate(
    plot = list(
      data |>
        ggplot(aes(TIMESTAMP, name, fill = name, alpha = coverage)) +
        geom_tile(height = 0.25, show.legend = FALSE) +
        labs(x = NULL, y = NULL, title = year)
    )
  ) |>
  pull(plot)

alt_cleaned <- alt_data |>
  # Clean vars with their QC flags
  map(
    \(d) map2(
      gf_vars, glue("QC_{gf_vars}"),
      \(x, y) transmute(d, "{x}" := if_else(.data[[y]] == 2, NA, .data[[x]]))
    )
  ) |>
  map(bind_cols) |>
  map(\(x) mutate(x, across(where(is.logical), as.numeric))) |>
  # Add level suffix to names
  map(\(x) rename_with(x, \(y) paste0(y, "_1"))) |>
  # Add back timestamp
  imap(\(x, idx) mutate(x, TIMESTAMP = alt_data[[idx]]$TIMESTAMP, .before = 1))

# Select focal ERA vars and add levels to single-level vars
era_cleaned <- era |>
  select(TIMESTAMP, starts_with(gf_vars)) |> 
  rename_with(\(x) paste0(x, "_1"), c(-TIMESTAMP, -matches("_\\d"))) |>
  # Add intermediate levels
  mutate(
    SWC_1.5 = (SWC_1 + SWC_2)/2, 
    SWC_2.5 = (SWC_2 + SWC_3)/2, 
    SWC_3.5 = (SWC_3 + SWC_4)/2,
    TS_1.5 = (TS_1 + TS_2)/2, 
    TS_2.5 = (TS_2 + TS_3)/2, 
    TS_3.5 = (TS_3 + TS_4)/2
  )

# Gather all alternate data
alt <- c(alt_cleaned, list(ERA = era_cleaned))

## De-bias alternate data ------------------------------------------------------

# TODO: separate loops for alt sites & ERA??
# - for some vars, alt sites are essentially drop-in replacements

# win_ref <- quarter(data$TIMESTAMP, "year.quarter")
win_ref <- factor(quarter(data$TIMESTAMP, "year.quarter"))
windows <- unique(win_ref)

debiased <- setNames(vector("list", length(gf_vars)), gf_vars)
fit_results <- setNames(vector("list", length(gf_vars)), gf_vars)

# setdiff(vars, soil_vars)
for (var in gf_vars) {
  
  # (var*window):source:level
  var_data <- cleaned[[var]]
  var_alt <- alt |>
    map(\(x) select(x, starts_with(var))) |>
    map(as.list) |>
    list_flatten()
  sources <- names(var_alt)

  debiased[[var]] <- vector("list", length(windows))
  fit_results[[var]] <- vector("list", length(windows))
  for (win in seq_along(windows)) {
    
    w_i <- which(as.integer(win_ref) == win)
    
    alt_fits <- setNames(vector("list", length(sources)), sources)
    alt_debiased <- setNames(vector("list", length(sources)), sources)
    for (src in sources) {
      var_src <- var_alt[[src]]
      
      # 1. Flexible window size if not enough good data
      win_src <- get_debias_window(var_src, var_data, w_i, na_frac)
      # Subset data within window
      var_win <- var_data[win_src]
      alt_win <- var_src[win_src]
      
      # 2. Remove phase difference
      lag_src <- get_best_lag(alt_win, var_win)
      alt_lag <- apply_lag(alt_win, lag_src)
      # Fill leading/trailing NAs (introduced by phase correction)
      alt_lag <- zoo::na.approx(alt_lag, rule = 2, maxgap = abs(lag_src))
      
      # 3. Correct bias - type of fit depends on var
      if (!any(complete.cases(var_win, alt_lag))) {
        # First check if fit exists
        fit <- NA
      } else if (var %in% cfg$biomet_gapfill$lm0_vars) {
        # Linear regression - forced through (0, 0)
        fit <- lm(var_win ~ 0 + alt_lag)
      } else if (var == "RH") {
        # Slightly different fit method for ERA vs. alt source
        if (str_starts(src, "ERA")) {
          # Regular linear regression
          fit <- lm(var_win ~ alt_lag)
        } else {
          # Linear regression - forced through (100, 100)
          fit <- lm(I(-var_win + 100) ~ 0 + I(-alt_lag + 100))
        }
      } else if (var == "WD") {
        # Offset that maximizes covariance
        fit <- cov_wd(alt_lag, var_win)
      } else {
        # Regular linear regression
        fit <- lm(var_win ~ alt_lag)
      }
      
      # Quality of fit and debiasing parameters
      fit_stats <- get_fit_stats(alt_lag, var_win, fit)
      fit_stats$lag <- lag_src
      fit_stats$n_fit <- length(win_src)
      fit_stats$n_win <- length(w_i)
      
      # Back-calculate intercept for RH fit (non-ERA alternates)
      if (var == "RH" & !str_starts(src, "ERA")) {
        fit_stats$int <- 100 * (1 - fit_stats$slope)
      }
      
      # Simple ratio used for precip slope
      if (var == "P_RAIN") {
        cc <- complete.cases(var_win, alt_lag)
        fit_stats$slope <- sum(var_win[cc])/sum(alt_lag[cc])
        fit_stats$int <- 0
      }
      
      # Debias phase-corrected alternate data
      alt_corr <- apply_lag(var_src, lag_src) * fit_stats$slope + fit_stats$int
      alt_corr <- alt_corr[w_i] # subset window
      # Fill leading/trailing NAs (introduced by phase correction)
      fit_pred <- zoo::na.approx(alt_corr, rule = 2, maxgap = abs(lag_src))
      
      # Upper limit on RH (ERA only - fit not forced through 100, 100)
      # TODO do the same with VPD??
      if (var == "RH" & str_starts(src, "ERA")) {
        fit_pred <- pmin(fit_pred, 100)
      }
      
      # Simplify WD degrees
      if (var == "WD") {
        fit_pred <- fit_pred %% 360
      }
      
      # Determine how many gaps were filled
      gap <- is.na(var_data[w_i])
      fit_stats$n_gap <- sum(gap)
      fit_stats$n_filled <- sum(gap & !is.na(fit_pred))
      fit_stats$window <- windows[win]
      
      alt_fits[[src]] <- fit_stats
      alt_debiased[[src]] <- fit_pred
      
    }
    
    # Compare alternate fits, determine fill order for window
    src_order <- order(map_dbl(alt_fits, "r2"), decreasing = TRUE)
    # # Prefer alternates that fill all gaps
    # all_filled <- map_lgl(alt_fits, \(x) x$n_filled == x$n_gap)
    # # Select alternate with highest R^2 in window
    # best_src <- which.max(map_dbl(alt_fits, "r2") * all_filled)
    # debiased[[var]][[win]] <- alt_debiased[[best_src]]
    # alt_fits[[best_src]]$n_used <- length(w_i)
    
    # Add predictions in order of fit quality
    db <- rep(NA, length(w_i))
    n_left <- length(w_i)
    for (i in src_order) {
      if (!anyNA(db)) break
      db <- coalesce(db, alt_debiased[[i]])
      alt_fits[[i]]$n_used <- n_left - sum(is.na(db))
      n_left <- sum(is.na(db)) # update # of NAs left to fill
    }
    debiased[[var]][[win]] <- db
    
    # Save results of all fits (for diagnostics & general reference)
    fit_results[[var]][[win]] <- alt_fits
  }
}

fit_results_df <- fit_results |>
  map_depth(2, \(x) bind_rows(x, .id = "name")) |> 
  map(\(x) bind_rows(x, .id = "name2")) |> 
  bind_rows(.id = "var") |>
  mutate(name = coalesce(name, name2)) |>
  select(-name2) |>
  mutate(name = str_remove(name, paste0(var, "_"))) |>
  relocate(window, .after = 1) |>
  separate_wider_delim(name, "_", names = c("src", "lvl"))

# fit_results_df |>
#   filter(n_used > 0)

debiased_df <- debiased |> 
  map(list_c) |>
  bind_cols() |>
  rename_with(\(x) paste0(x, "_ALT"))

cleaned |>
  mutate(G = zoo::na.approx(G, na.rm = FALSE, maxgap = max_interp)) |>
  ggplot(aes(TIMESTAMP, coalesce(G, debiased_df$G_ALT))) +
  geom_line(color = "coral") +
  geom_line(aes(y = G))
cleaned |>
  mutate(TS = zoo::na.approx(TS, na.rm = FALSE, maxgap = max_interp)) |>
  ggplot(aes(TIMESTAMP, coalesce(TS, debiased_df$TS_ALT))) +
  geom_line(color = "coral") +
  geom_line(aes(y = TS))
cleaned |>
  mutate(SWC = zoo::na.approx(SWC, na.rm = FALSE, maxgap = max_interp)) |>
  ggplot(aes(TIMESTAMP, coalesce(SWC, debiased_df$SWC_ALT))) +
  geom_line(color = "coral") +
  geom_line(aes(y = SWC))


## Time series-based gap-filling methods ---------------------------------------

### Linear interpolation ----

interp_filled <- vector("list", length(gf_vars))
names(interp_filled) <- paste0(gf_vars, "_LI")
for (i in seq_along(gf_vars)) {
  interp_filled[[i]] <- zoo::na.approx(
    cleaned[[gf_vars[i]]], na.rm = FALSE, maxgap = max_interp
  )
}

### Kalman smoothing on state space representation of ARIMA model ----

# For highly autocorrelated site-specific soil vars 

# kalman_filled <- soil_vars |>
#   map(
#     \(x) imputeTS::na_kalman(
#       cleaned[[x]], model = "auto.arima", maxgap = max_kalman, stepwise = FALSE
#     )
#   ) |>
#   set_names(glue("{soil_vars}_KS")) |>
#   bind_cols()

### Mean diurnal course & look up table ----

mds_vars <- vars |> 
  keep(\(x) "MDS" %in% x$gapfill) |>
  names()

cleaned_reddyproc <- cleaned |>
  mutate(TIMESTAMP = TIMESTAMP + minutes(15)) |>
  as.data.frame()

mds <- sEddyProc$new(site, cleaned_reddyproc, mds_vars, "TIMESTAMP")

for (var in mds_vars) {
  mds$sMDSGapFill(
    var, 
    QFVar = "none", 
    QFValue = NA, 
    V1 = "SW_IN", 
    V2 = "VPD", 
    V3 = "TA",
    FillAll = FALSE, 
    isVerbose = FALSE
  )
}

# Set filled vars to NA if poor quality (FQC = 2 or 3)
mds_results <- mds$sExportResults()
mds_filled <- vector("list", length(mds_vars))
names(mds_filled) <- paste0(mds_vars, "_MDS")
for (i in seq_along(mds_vars)) {
  mds_filled[[i]] <- if_else(
    mds_results[[paste0(mds_vars[i], "_fqc")]] > 1, 
    NA, 
    mds_results[[paste0(mds_vars[i], "_f")]]
  )
}

## Fill gaps -------------------------------------------------------------------

### Gather de-biased & filled vars ----

filled <- bind_cols(
  cleaned, debiased_df, bind_cols(interp_filled), bind_cols(mds_filled)
)
for (var in gf_vars) {
  fill_seq <- c(var, paste0(var, "_", vars[[var]]$gapfill))
  fill_data <- as.list(filled[, fill_seq])
  var_filled <- coalesce(!!!fill_data)
  filled[[paste0(var, "_F")]] <- var_filled
}
filled <- select(filled, ends_with("_F"))

nrow(filled) == nrow(drop_na(filled)) # should be TRUE
map_lgl(filled, anyNA) # should all be FALSE

filled <- bind_cols(data, filled)

# Update nighttime indicator (REddyProc definition)
filled$NIGHT <- filled$SW_IN_F <= 10 & filled$SW_IN_POT == 0

# Update LI-7500 QC flags with gap-filled precipitation
data$QC_PRECIP_F = if_else(filled$P_RAIN_F >= precip_thr, 2, 0)
filled$QC_LE <- combn_QC(as.data.frame(data), c("QC_LE", "QC_PRECIP_F"))
filled$QC_FC <- combn_QC(as.data.frame(data), c("QC_FC", "QC_PRECIP_F"))

# Save PDF of time series plots showing gap-filled sections
report_file <- glue("metgf_plots_{site}.pdf")
pdf(file.path(report_dir, report_file))
coverage_plot
map(
  gf_vars,
  \(x) cleaned |>
    ggplot(aes(TIMESTAMP)) +
    geom_line(
      aes(y = filled[[glue("{x}_F")]], color = "GF"), 
      show.legend = FALSE, linewidth = 0.3
    ) +
    geom_line(aes(y = .data[[x]]), linewidth = 0.3, na.rm = TRUE) +
    labs(y = x, x = NULL)
)
dev.off()

## Write gap-filled dataset with documentation ---------------------------------

# Separately by year
for (i in seq_along(years)) {
  gf_file <- glue("eddypro_{site}-{years[i]}_fluxnet_adv_metgf.csv")
  write_csv(
    filter(filled, year(TIMESTAMP) == years[i]), 
    file.path(output_dirs[i], gf_file)
  )
}

# Documentation
# [goes here]

