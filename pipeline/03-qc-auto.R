# Automatic quality control ----------------------------------------------------

# Purpose: 

# Load the required packages
library("remotes")
if (!require("openeddy")) {
  install_github("lsigut/openeddy")
}
library("yaml")
library("glue")
library("tidyverse")

# Load custom functions
source("R/pipeline-funs.R")
source("R/flag-spikes.R")

## Initialize script settings & documentation ----------------------------------

data_dir <- file.path("data", site) # path to directory containing all data
output_dir <- file.path("output", site, year) # where to write output
report_dir <- file.path("reports", site) # where to save reports

# Load metadata file
md <- read_yaml(file.path(data_dir, "metadata.yml"))

# Load general configuration
cfg <- read_yaml("config/config.yml")

# Load variable configuration
vars <- read_yaml("config/var-config.yml")

# Set configuration options
mvc_dist_thr <- cfg$qc_auto$mvc_dist_thr
rad_dist_thr <- cfg$qc_auto$rad_dist_thr
night_thr <- cfg$qc_auto$night_thr
agc_min <- cfg$qc_auto$agc_min
rssi_min <- cfg$qc_auto$rssi_min
scf_max <- cfg$qc_auto$scf_max
rain_thr <- cfg$qc_auto$rain_thr

## Read corrected data file ----------------------------------------------------

data_file <- glue("eddypro_{site}-{year}_fluxnet_adv_corrected.csv")

data <- read_csv(file.path(output_dir, data_file))

## Create diagnostic plots for manual inspection -------------------------------

# Select subset of variables for diagnostics
dep_vars <- vars[cfg$dep_vars]

flux_vars <- vars[cfg$flux_vars]

# Save PDF of time series plots showing dependencies

dep_report_file <- glue("qc_dependencies_{site}-{year}.pdf")
pdf(file.path(report_dir, dep_report_file), width = 11, height = 8.5)

for (var in names(dep_vars)) {
  data |> 
    rename(timestamp = TIMESTAMP) |> 
    mutate("{var}" := structure(.data[[var]], units = vars[[var]]$units)) |>
    as.data.frame() |> 
    plot_precheck(var)
}

dev.off()

## Perform QC checks on dependencies -------------------------------------------

cleaned <- data

### Physical limits ----

# Dependencies are cleaned directly

# Get dependencies with defined plausible limits
dep_lims <- map(dep_vars, "phys_lims")

# Remove values outside limits
for (var in names(dep_lims)) {
  lim <- dep_lims[[var]]
  cleaned[[var]] <- na_if_outside(cleaned[[var]], lim[1], lim[2])
}

### Spikes ----

# [spike detection for dependencies?]

### Combine dependencies by instrument ----

sa_vars <- cfg$qc_auto$sa_vars
# Flag if NA introduced during dependency cleaning
sa_dep <- apply(is.na(cleaned[sa_vars]) & !is.na(data[sa_vars]), 1, any)
# Restore original NAs
sa_dep[apply(is.na(data[sa_vars]), 1, any)] <- NA

li7500_vars <- cfg$qc_auto$li7500_vars
# Flag if NA introduced during dependency cleaning
li7500_dep <- apply(
  is.na(cleaned[li7500_vars]) & !is.na(data[li7500_vars]), 1, any
)
# Restore original NAs
li7500_dep[apply(is.na(data[li7500_vars]), 1, any)] <- NA

li7700_vars <- cfg$qc_auto$li7700_vars
# Flag if NA introduced during dependency cleaning
li7700_dep <- apply(
  is.na(cleaned[li7700_vars]) & !is.na(data[li7700_vars]), 1, any
)
# Restore original NAs
li7700_dep[apply(is.na(data[li7700_vars]), 1, any)] <- NA

## Perform QC checks on biomet data --------------------------------------------

met_qc <- data[, 0]

### Dependencies ----

wind_dep <- is.na(cleaned$U_SIGMA + cleaned$V_SIGMA) & !is.na(data$U_SIGMA)

met_qc$QC_WS_DEP <- if_else(sa_dep | wind_dep, 2, 0)
met_qc$QC_WD_DEP <- met_qc$QC_WS_DEP
met_qc$QC_USTAR_DEP <- met_qc$QC_WS_DEP

### Physical limits ----

# Get vars with defined plausible limits (excluding dependencies)
met_lims <- vars |> 
  keep(\(x) str_starts(x$type, "met")) |> 
  map("phys_lims") |> 
  compact()

# Flag values outside limits
for (var in names(met_lims)) {
  lims <- met_lims[[var]]
  flag <- apply_thr(as.numeric(data[[var]]), lims, flag = "outside")
  met_qc[[glue("QC_{var}_LIM")]] <- flag
}

### Spikes ----

spike_vars <- cfg$qc_auto$spike_vars

for (var in spike_vars) {
  flag <- flag_spikes(data[[var]], lim_flag = met_qc[[glue("QC_{var}_LIM")]])
  met_qc[[glue("QC_{var}_SPIKE")]] <- flag
}

### Multivariate comparison ----

# TODO: save plots to a file
mvc_plots <- list()

# SW_IN vs PPFD_IN
if (any(complete.cases(data$PPFD_IN, data$SW_IN))) {
  rad_odr <- odr_coef(data$PPFD_IN, data$SW_IN)
  rad_resid <- data$SW_IN - (data$PPFD_IN * rad_odr[1] + rad_odr[2])
  rad_dist <- abs(rad_resid/sqrt(1 + rad_odr[1]^2))
  rad_rse <- sqrt(mean(rad_dist^2, na.rm = TRUE))
  met_qc$QC_SW_IN_MVC <- if_else(
    rad_dist > mvc_dist_thr * rad_rse & rad_dist > rad_dist_thr, 2, 0
  )
} else {
  met_qc$QC_SW_IN_MVC <- rep(NA, nrow(data))
}
met_qc$QC_PPFD_IN_MVC <- met_qc$QC_SW_IN_MVC

# Save MVC plot
mvc_plots$rad <- data |>
  mutate(flag = factor(met_qc$QC_SW_IN_MVC)) |>
  ggplot(aes(PPFD_IN, SW_IN)) +
  geom_point(aes(color = flag), na.rm = TRUE) +
  geom_smooth(method = "lm", se = FALSE, color = "black", na.rm = TRUE)

# WS vs USTAR
temp_ws <- if_else(met_qc$QC_WS_DEP == 2, NA_real_, data$WS)
temp_ustar <- if_else(met_qc$QC_USTAR_DEP == 2, NA_real_, data$USTAR)
if (any(complete.cases(temp_ustar, temp_ws))) {
  wind_odr <- odr_coef(temp_ustar, temp_ws)
  wind_resid <- temp_ws - (temp_ustar * wind_odr[1] + wind_odr[2])
  wind_dist <- abs(wind_resid/sqrt(1 + wind_odr[1]^2))
  wind_rse <- sqrt(mean(wind_dist^2, na.rm = TRUE))
  met_qc$QC_WS_MVC <- if_else(wind_dist > mvc_dist_thr * wind_rse, 2, 0)
} else {
  met_qc$QC_WS_MVC <- rep(NA, nrow(data))
}
# NA flags are 0 since used cleaned vars for test
# met_qc$QC_WS_MVC <- replace_na(met_qc$QC_WS_MVC, 0)
met_qc$QC_USTAR_MVC <- met_qc$QC_WS_MVC

# Save MVC plot
mvc_plots$wind <- cleaned |>
  mutate(
    flag = factor(met_qc$QC_WS_MVC),
    WS = temp_ws,
    USTAR = temp_ustar
  ) |>
  ggplot(aes(USTAR, WS)) +
  geom_point(aes(color = flag), na.rm = TRUE) +
  geom_smooth(method = "lm", se = FALSE, color = "black", na.rm = TRUE)

# TA vs T_SONIC
if (any(complete.cases(cleaned$T_SONIC, data$TA))) {
  temp_odr <- odr_coef(cleaned$T_SONIC, data$TA)
  temp_resid <- data$TA - (cleaned$T_SONIC * temp_odr[1] + temp_odr[2])
  temp_dist <- abs(temp_resid/sqrt(1 + temp_odr[1]^2))
  temp_rse <- sqrt(mean(temp_dist^2, na.rm = TRUE))
  met_qc$QC_TA_MVC <- if_else(temp_dist > mvc_dist_thr * temp_rse, 2, 0)
  # met_qc$QC_T_SONIC_MVC <- met_qc$QC_TA_MVC
} else {
  met_qc$QC_TA_MVC <- rep(NA, nrow(data))
}
# NA flags are 0 since used cleaned vars for test
# met_qc$QC_TA_MVC <- replace_na(met_qc$QC_TA_MVC, 0)

# Save MVC plot
mvc_plots$temp <- cleaned |>
  mutate(flag = factor(met_qc$QC_TA_MVC)) |>
  ggplot(aes(T_SONIC, TA)) +
  geom_point(aes(color = flag), na.rm = TRUE) +
  geom_smooth(method = "lm", se = FALSE, color = "black", na.rm = TRUE)


## Perform QC checks on flux data ----------------------------------------------

flux_qc <- data[, 0]

### Dependencies ----

flux_qc$QC_H_DEP <- if_else(sa_dep, 2, 0)
flux_qc$QC_LE_DEP <- if_else(sa_dep | li7500_dep, 2, 0)
flux_qc$QC_FC_DEP <- if_else(sa_dep | li7500_dep, 2, 0)
flux_qc$QC_FCH4_DEP <- if_else(sa_dep | li7700_dep, 2, 0)

### SSITC flags ----

flux_qc$QC_H_SSITC <- data$H_SSITC_TEST
flux_qc$QC_LE_SSITC <- data$LE_SSITC_TEST
flux_qc$QC_FC_SSITC <- data$FC_SSITC_TEST
flux_qc$QC_FCH4_SSITC <- data$FCH4_SSITC_TEST

### Plausible time lag ----

flux_qc$QC_LE_TLAG <- if_else(data$H2O_TLAG_ACTUAL != data$H2O_TLAG_USED, 2, 0)
flux_qc$QC_FC_TLAG <- if_else(data$CO2_TLAG_ACTUAL != data$CO2_TLAG_USED, 2, 0)
flux_qc$QC_FCH4_TLAG <- if_else(data$CH4_TLAG_ACTUAL != data$CH4_TLAG_USED, 2, 0)

### Spectral correction factor ----

scf_thr <- rep(scf_max, 2)
flux_qc$QC_H_SCF <- apply_thr(data$H_SCF, scf_thr, flag = "higher")
flux_qc$QC_LE_SCF <- apply_thr(data$LE_SCF, scf_thr, flag = "higher")
flux_qc$QC_FC_SCF <- apply_thr(data$FC_SCF, scf_thr, flag = "higher")
flux_qc$QC_FCH4_SCF <- apply_thr(data$FCH4_SCF, scf_thr, flag = "higher")

### Sensor signal strength ----

# LI-7500: LE, FC
agc_flag <- apply_thr(
  data$INST_LI7500_AGC_OR_RSSI, rep(agc_min, 2), flag = "lower"
)
flux_qc$QC_LE_AGC <- agc_flag
flux_qc$QC_FC_AGC <- agc_flag

# LI-7700: FCH4
flux_qc$QC_FCH4_RSSI <- apply_thr(
  data$CUSTOM_RSSI_77_MEAN, rep(rssi_min, 2), flag = "lower"
)

### Precipitation ----

precip_flag <- apply_thr(
  as.numeric(data$P_RAIN), rep(rain_thr, 2), flag = "higher"
)
flux_qc$QC_LE_PRECIP <- precip_flag
flux_qc$QC_FC_PRECIP <- precip_flag

### Low frequency spikes ----

# Make preliminary flags to apply before despiking
flux_qc$QC_H_PRELIM <- combn_QC(
  as.data.frame(flux_qc), 
  qc_names = str_subset(names(flux_qc), "_H_"), 
  no_messages = TRUE
)
flux_qc$QC_LE_PRELIM <- combn_QC(
  as.data.frame(flux_qc), 
  qc_names = str_subset(names(flux_qc), "_LE_"),
  na.as_0_pattern = "PRECIP$", 
  no_messages = TRUE
)
flux_qc$QC_FC_PRELIM <- combn_QC(
  as.data.frame(flux_qc), 
  qc_names = str_subset(names(flux_qc), "_FC_"),
  na.as_0_pattern = "PRECIP$", 
  no_messages = TRUE
)

# 'despikeLF' expects specific naming
flux_qc$timestamp <- data$TIMESTAMP
flux_qc$GR <- data$SW_IN_POT

flux_qc$QC_H_SPIKESLF <- despikeLF(
  as.data.frame(bind_cols(data, flux_qc)), 
  var = "H", 
  qc_flag = "QC_H_PRELIM", 
  var_thr = vars$H$phys_lims, 
  light = "GR", 
  night_thr = night_thr
)

flux_qc$QC_LE_SPIKESLF <- despikeLF(
  as.data.frame(bind_cols(data, flux_qc)), 
  var = "LE", 
  qc_flag = "QC_LE_PRELIM", 
  var_thr = vars$LE$phys_lims, 
  light = "GR", 
  night_thr = night_thr
)

flux_qc$QC_FC_SPIKESLF <- despikeLF(
  as.data.frame(bind_cols(data, flux_qc)), 
  var = "FC", 
  qc_flag = "QC_FC_PRELIM", 
  var_thr = vars$FC$phys_lims, 
  light = "GR", 
  night_thr = night_thr
)

# FCH4 is checked for range but not spikes
fch4_lims <- cfg$qc_auto$fch4_lims[[site]]
flux_qc$QC_FCH4_LIM <- apply_thr(data$FCH4, fch4_lims, flag = "outside")

flux_qc <- select(flux_qc, -timestamp, -GR, -ends_with("_PRELIM"))

### Footprint (WIP) ----

# source("R/footprint-model.R")
# 
# # Subset variables for fetch calculations
# data_fetch <- cleaned |>
#   select(
#     ws = WS, ustar = USTAR, mo_length = MO_LENGTH, z = BADM_INST_HEIGHT_SA,
#     zd = DISPLACEMENT_HEIGHT, blh = BLH_ERA
#   ) |>
#   mutate(pct = 0.85, .before = 1)
# 
# # Calculate fetch
# fetch <- data_fetch |>
#   pmap(calc_fetch_FFP, .progress = TRUE) |>
#   bind_rows() |>
#   mutate(TIMESTAMP = cleaned$TIMESTAMP, .before = 1)


### Summary plots ----

# [summary plots go here]

## Combine flags ---------------------------------------------------------------

met_vars <- unique(str_match(names(met_qc), "QC_(.*)_[^_]*$")[, 2])
met_flags <- met_vars |> 
  set_names(glue("QC_{met_vars}")) |> 
  map(\(x) str_subset(names(met_qc), glue("_{x}_"))) |> 
  map(
    \(nm) combn_QC(
      x = as.data.frame(met_qc), 
      qc_names = nm, 
      additive = FALSE, 
      # na.as = NA, 
      # For MVC flags NA = 0 since used cleaned vars for test
      na.as_0_pattern = "MVC$"
    )
  ) |> 
  bind_cols()

flux_vars <- c("H", "LE", "FC", "FCH4")
flux_flags <- flux_vars |> 
  set_names(glue("QC_{flux_vars}")) |> 
  map(\(x) str_subset(names(flux_qc), glue("_{x}_"))) |> 
  map(
    \(x) combn_QC(
      as.data.frame(flux_qc), x, additive = FALSE, 
      na.as_0_pattern = "SPIKESLF$|PRECIP$"
    )
  ) |> 
  bind_cols()

cleaned <- bind_cols(cleaned, met_flags, flux_flags)

met_qc_summary <- met_qc |> 
  summarize(across(everything(), \(x) sum(x > 1, na.rm = TRUE))) |>
  pivot_longer(
    everything(), names_to = c("var", "test"), 
    names_pattern = "QC_(.+_*.*)_(.+)"
  ) |>
  pivot_wider(names_from = test, values_from = value)

flux_qc_summary <- flux_qc |> 
  summarize(across(everything(), \(x) sum(x > 1, na.rm = TRUE))) |>
  pivot_longer(
    everything(), names_to = c("flux", "test"), names_pattern = "QC_(.+)_(.+)"
  ) |>
  pivot_wider(names_from = test, values_from = value)

# summary_QC(
#   as.data.frame(flux_qc),
#   c("QC_FC_DEP", "QC_FC_SSITC", "QC_FC_TLAG", "QC_FC_SCF", "QC_FC_AGC", "QC_FC_PRECIP", "QC_FC_SPIKESLF"),
#   na.as_0_pattern = "SPIKESLF$|PRECIP$",
#   cumul = TRUE, plot = TRUE
# )

# TODO: save these summaries


## Summary plots ---------------------------------------------------------------

### Biomet variables ----

biomet_report_file <- glue("qc_biomet_{site}-{year}.pdf")
pdf(file.path(report_dir, biomet_report_file), width = 11, height = 8.5)

for (var in met_vars) {
  # Workaround for plot_eddy error when no valid data
  if (all(is.na(cleaned[[var]]))) next
  
  cleaned |> 
    rename(timestamp = TIMESTAMP) |> 
    mutate("{var}" := structure(.data[[var]], units = vars[[var]]$units)) |> 
    as.data.frame() |> 
    plot_eddy(
      var, 
      qc_flag = glue("QC_{var}"),
      test = glue("QC_{var}"),
      skip = "weekly",
      ylim = vars[[var]]$phys_lims,
      document = TRUE
    )
}

dev.off()

### Fluxes ----

flux_flag_vars <- flux_vars |> 
  map(\(x) str_subset(names(flux_qc), glue("_{x}_"))) |> 
  set_names(flux_vars)

flux_report_file <- glue("qc_fluxes_{site}-{year}.pdf")
pdf(file.path(report_dir, flux_report_file), width = 11, height = 8.5)

# QC flag summaries
for (var in flux_vars) {
  gridExtra::grid.arrange(
    grobs = list(
      summary_QC(
        as.data.frame(flux_qc),
        flux_flag_vars[[var]],
        na.as_0_pattern = "SPIKESLF$|PRECIP$",
        cumul = FALSE, plot = TRUE
      ),
      summary_QC(
        as.data.frame(flux_qc),
        flux_flag_vars[[var]],
        na.as_0_pattern = "SPIKESLF$|PRECIP$",
        cumul = TRUE, plot = TRUE
      )
    ),
    nrow = 1
  )
}

# Diagnostic time series plots
for (var in flux_vars) {
  data |> 
    rename(
      timestamp = TIMESTAMP, PAR = PPFD_IN, GR = SW_IN, Rn = NETRAD,
      Tair = TA, Tsoil = TS
    ) |> 
    mutate(
      PAR = structure(PAR, units = vars[["PPFD_IN"]]$units),
      GR = structure(GR, units = vars[["SW_IN"]]$units),
      VPD = structure(VPD, units = vars[["VPD"]]$units),
      Rn = structure(Rn, units = vars[["NETRAD"]]$units),
      Tair = structure(Tair, units = vars[["TA"]]$units),
      Tsoil = structure(Tsoil, units = vars[["TS"]]$units),
      "{var}" := structure(.data[[var]], units = vars[[var]]$units)
    ) |>
    as.data.frame() |> 
    plot_eddy(
      var, 
      qc_flag = glue("{var}_SSITC_TEST"), 
      test = glue("{var}_SSITC_TEST"),
      ylim = vars[[var]]$phys_lims,
      document = TRUE
    )
}

dev.off()

## Write quality-flagged dataset with documentation ----------------------------

qc_file <- glue("eddypro_{site}-{year}_fluxnet_adv_qaqc.csv")
write_csv(cleaned, file.path(output_dir, qc_file))

# Documentation
# [goes here]

