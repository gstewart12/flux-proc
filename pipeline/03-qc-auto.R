# Automatic quality control ----------------------------------------------------

# Purpose: 

# TODO: store parameters in external file

# Load the required packages
library("remotes")
if (!require("openeddy")) {
  install_github("lsigut/openeddy")
}
library("yaml")
library("glue")
library("tidyverse")

source("source/flag-spikes.R")

## Initialize script settings & documentation ----------------------------------

site <- "JLR" # three-letter site code
year <- 2023 # a single year to process
output_dir <- file.path("output", site, year) # where to write output

# Load metadata file
md <- read_yaml(file.path("data", site, "metadata.yml"))

# Load variable info file
vars <- read_yaml("var-info.yml")

# Set configuration options
mvc_dist_thr <- 5
rad_dist_thr <- 50
night_thr <- 20 # change to 50?
agc_min <- 70 # min LI-7500 (CO2/H2O) signal strength
rssi_min <- 10 # min LI-7700 (CH4) signal strength
scf_max <- 2 # max spectral correction factor
rain_thr <- 0 # max rainfall for data depending on LI-7500

## Read corrected data file ----------------------------------------------------

data_file <- glue("eddypro_{site}-{year}_fluxnet_adv_corrected.csv")

data <- read_csv(file.path(output_dir, data_file))

## Perform QC checks on dependencies -------------------------------------------

cleaned <- data

### Plausible limits ----

na_if_outside <- function(x, left, right) {
  if_else(between(x, left, right), x, NA_real_)
}

# Dependencies are cleaned directly

# TODO: load these from external file

cleaned <- cleaned |>
  mutate(
    T_SONIC = na_if_outside(T_SONIC, -20, 50),
    T_SONIC_SIGMA = na_if_outside(T_SONIC_SIGMA, 0, 2),
    U_UNROT = na_if_outside(U_UNROT, -15, 15),
    U_SIGMA = na_if_outside(U_SIGMA, 0, 5),
    V_UNROT = na_if_outside(V_UNROT, -15, 15),
    V_SIGMA = na_if_outside(V_SIGMA, 0, 3),
    W_UNROT = na_if_outside(W_UNROT, -0.7, 0.7),
    W_SIGMA = na_if_outside(W_SIGMA, 0, 1.5),
    ROT_PITCH = na_if_outside(ROT_PITCH, -30, 30), # second axis rotation angle
    H2O = na_if_outside(H2O, 0, 45),
    H2O_SIGMA = na_if_outside(H2O_SIGMA, 0, 100),
    CO2 = na_if_outside(CO2, 250, 1000),
    CO2_SIGMA = na_if_outside(CO2_SIGMA, 0, 1.5),
    CH4 = na_if_outside(CH4, 0, 15),
    CH4_SIGMA = na_if_outside(CH4_SIGMA, 0, 0.025)
  )

### Spikes ----

# [spike detection for dependencies?]

### Combine dependencies by instrument ----

sa_vars <- c(
  "U_UNROT", "U_SIGMA", "V_UNROT", "V_SIGMA", "W_UNROT", "W_SIGMA", 
  "T_SONIC", "T_SONIC_SIGMA", "ROT_PITCH"
)
sa_dep <- sa_vars |>
  map(\(x) is.na(cleaned[[x]]) & !is.na(data[[x]])) |>
  reduce(`|`)

li7500_vars <- c("H2O", "H2O_SIGMA", "CO2", "CO2_SIGMA")
li7500_dep <- li7500_vars |>
  map(\(x) is.na(cleaned[[x]]) & !is.na(data[[x]])) |>
  reduce(`|`)

li7700_vars <- c("CH4", "CH4_SIGMA")
li7700_dep <- li7700_vars |>
  map(\(x) is.na(cleaned[[x]]) & !is.na(data[[x]])) |>
  reduce(`|`)

## Perform QC checks on biomet data --------------------------------------------

# Met vars: TA, PA, RH, VPD, WS, WD, USTAR 
# Biomet vars: PPFD_IN, SW_IN, SW_OUT, LW_IN, LW_OUT, NETRAD, P_RAIN, G, TS, SWC

met_qc <- data[, 0]

### Dependencies ----

wind_dep <- is.na(cleaned$U_SIGMA + cleaned$V_SIGMA) & !is.na(data$U_SIGMA)

met_qc$QC_WS_DEP <- if_else(sa_dep | wind_dep, 2, 0)
met_qc$QC_WD_DEP <- met_qc$QC_WS_DEP
met_qc$QC_USTAR_DEP <- met_qc$QC_WS_DEP

### Plausible limits ----

# vars |> map("phys_lims") |> compact()

met_qc <- mutate(
  met_qc,
  QC_TA_LIM = if_else(!between(data$TA, -20, 50), 2, 0),
  QC_PA_LIM = if_else(!between(data$PA, 95, 105), 2, 0),
  QC_RH_LIM = if_else(!between(data$RH, 0, 100), 2, 0),
  QC_VPD_LIM = if_else(!between(data$VPD, 0, 45), 2, 0),
  QC_WS_LIM = if_else(!between(data$WS, 0, 40), 2, 0),
  QC_WD_LIM = if_else(!between(data$WD, 0, 360), 2, 0),
  QC_USTAR_LIM = if_else(!between(data$USTAR, 0, 3), 2, 0),
  QC_PPFD_IN_LIM = if_else(!between(data$PPFD_IN, 0, 2200), 2, 0),
  QC_SW_IN_LIM = if_else(!between(data$SW_IN, 0, 1200), 2, 0),
  QC_SW_OUT_LIM = if_else(!between(data$SW_OUT, -10, 500), 2, 0),
  QC_LW_IN_LIM = if_else(!between(data$LW_IN, 100, 600), 2, 0),
  QC_LW_OUT_LIM = if_else(!between(data$LW_OUT, 100, 750), 2, 0),
  QC_NETRAD_LIM = if_else(!between(data$NETRAD, -200, 900), 2, 0),
  QC_G_LIM = if_else(!between(data$G, -100, 250), 2, 0),
  QC_TS_LIM = if_else(!between(data$TS, -5, 35), 2, 0),
  QC_SWC_LIM = if_else(!between(data$SWC, 0, 100), 2, 0),
  QC_P_RAIN_LIM = if_else(!between(data$P_RAIN, 0, 50), 2, 0)
)

### Spikes ----

met_qc <- mutate(
  met_qc,
  QC_TA_SPIKE = flag_spikes(data$TA, lim_flag = QC_TA_LIM),
  QC_PA_SPIKE = flag_spikes(data$PA, lim_flag = QC_PA_LIM),
  QC_RH_SPIKE = flag_spikes(data$RH, lim_flag = QC_RH_LIM),
  QC_VPD_SPIKE = flag_spikes(data$VPD, lim_flag = QC_VPD_LIM),
  QC_LW_IN_SPIKE = flag_spikes(data$LW_IN, lim_flag = QC_LW_IN_LIM),
  QC_LW_OUT_SPIKE = flag_spikes(data$LW_OUT, lim_flag = QC_LW_OUT_LIM),
  QC_G_SPIKE = flag_spikes(data$G, lim_flag = QC_G_LIM),
  # Not checking TS spikes - too sensitive
  # QC_TS_SPIKE = flag_spikes(data$TS, lim_flag = QC_TS_LIM),
  QC_SWC_SPIKE = flag_spikes(data$SWC, lim_flag = QC_SWC_LIM)
)

### Multivariate comparison ----

odr_coef <- function(x, y) {
  
  # Orthogonal distance regression
  df <- data.frame(x = x, y = y)
  r <- prcomp(~ x + y, data = df, na.action = na.exclude)
  slope <- r$rotation[2, 1] / r$rotation[1, 1]
  int <- r$center[2] - slope * r$center[1]
  
  c(slope, int)
}

flag_mvc <- function(x, y, dist_thr, abs_dist_thr = 0) {
  odr <- odr_coef(x, y)
  resid <- y - (x * odr[1] + odr[2])
  dist <- abs(resid/sqrt(1 + odr[1]^2))
  rse <- sqrt(mean(dist^2, na.rm = TRUE))
  dplyr::if_else(dist > dist_thr * rse & dist > abs_dist_thr, 2, 0)
}

# SW_IN vs PPFD_IN
rad_odr <- odr_coef(data$PPFD_IN, data$SW_IN)
rad_resid <- data$SW_IN - (data$PPFD_IN * rad_odr[1] + rad_odr[2])
rad_dist <- abs(rad_resid/sqrt(1 + rad_odr[1]^2))
rad_rse <- sqrt(mean(rad_dist^2, na.rm = TRUE))
met_qc$QC_SW_IN_MVC <- if_else(
  rad_dist > mvc_dist_thr * rad_rse & rad_dist > rad_dist_thr, 2, 0
)
met_qc$QC_PPFD_IN_MVC <- met_qc$QC_SW_IN_MVC

# data |> 
#   mutate(flag = factor(met_qc$QC_SW_IN_MVC)) |> 
#   ggplot(aes(PPFD_IN, SW_IN, color = flag)) + 
#   geom_point()

# WS vs USTAR
temp_ws <- if_else(met_qc$QC_WS_DEP == 2, NA_real_, data$WS)
temp_ustar <- if_else(met_qc$QC_USTAR_DEP == 2, NA_real_, data$USTAR)
wind_odr <- odr_coef(temp_ustar, temp_ws)
wind_resid <- temp_ws - (temp_ustar * wind_odr[1] + wind_odr[2])
wind_dist <- abs(wind_resid/sqrt(1 + wind_odr[1]^2))
wind_rse <- sqrt(mean(wind_dist^2, na.rm = TRUE))
met_qc$QC_WS_MVC <- if_else(wind_dist > mvc_dist_thr * wind_rse, 2, 0)
# NA flags are 0 since used cleaned vars for test
met_qc$QC_WS_MVC <- replace_na(met_qc$QC_WS_MVC, 0)
met_qc$QC_USTAR_MVC <- met_qc$QC_WS_MVC

# cleaned |> 
#   mutate(flag = factor(met_qc$QC_WS_MVC)) |> 
#   ggplot(aes(USTAR, WS, color = flag)) + 
#   geom_point()

# TA vs T_SONIC
temp_odr <- odr_coef(data$T_SONIC, data$TA)
temp_resid <- data$TA - (data$T_SONIC * temp_odr[1] + temp_odr[2])
temp_dist <- abs(temp_resid/sqrt(1 + temp_odr[1]^2))
# temp_rse <- sqrt(mean(temp_resid^2, na.rm = TRUE))
temp_rse <- sqrt(mean(temp_dist^2, na.rm = TRUE))
met_qc$QC_TA_MVC <- if_else(temp_dist > mvc_dist_thr * temp_rse, 2, 0)
# met_qc$QC_T_SONIC_MVC <- met_qc$QC_TA_MVC

# cleaned |> 
#   mutate(flag = factor(met_qc$QC_TA_MVC)) |> 
#   ggplot(aes(T_SONIC, TA, color = flag)) + 
#   geom_point()

### Summary plots ----

# [summary plots go here]

## Perform QC checks on flux data ----------------------------------------------

flux_qc <- data[, 0]

### Dependencies ----

flux_qc <- mutate(
  flux_qc,
  QC_H_DEP = if_else(sa_dep, 2, 0),
  QC_LE_DEP = if_else(sa_dep | li7500_dep, 2, 0),
  QC_FC_DEP = if_else(sa_dep | li7500_dep, 2, 0),
  QC_FCH4_DEP = if_else(sa_dep | li7700_dep, 2, 0)
)

### SSITC flags ----

flux_qc <- mutate(
  flux_qc,
  QC_H_SSITC = data$H_SSITC_TEST,
  QC_LE_SSITC = data$LE_SSITC_TEST,
  QC_FC_SSITC = data$FC_SSITC_TEST,
  QC_FCH4_SSITC = data$FCH4_SSITC_TEST
)

### Other ----

# Plausible time lag

flux_qc <- mutate(
  flux_qc,
  QC_LE_TLAG = if_else(data$H2O_TLAG_ACTUAL != data$H2O_TLAG_USED, 2, 0),
  QC_FC_TLAG = if_else(data$CO2_TLAG_ACTUAL != data$CO2_TLAG_USED, 2, 0),
  QC_FCH4_TLAG = if_else(data$CH4_TLAG_ACTUAL != data$CH4_TLAG_USED, 2, 0)
)

# Spectral correction factor

flux_qc <- mutate(
  flux_qc,
  QC_H_SCF = if_else(data$H_SCF > scf_max, 2, 0),
  QC_LE_SCF = if_else(data$LE_SCF > scf_max, 2, 0),
  QC_FC_SCF = if_else(data$FC_SCF > scf_max, 2, 0),
  QC_FCH4_SCF = if_else(data$FCH4_SCF > scf_max, 2, 0)
)

# Sensor signal strength

flux_qc <- mutate(
  flux_qc,
  QC_LE_AGC = if_else(data$INST_LI7500_AGC_OR_RSSI < agc_min, 2, 0),
  QC_FC_AGC = if_else(data$INST_LI7500_AGC_OR_RSSI < agc_min, 2, 0),
  QC_FCH4_RSSI = if_else(data$CUSTOM_RSSI_77_MEAN < rssi_min, 2, 0)
)

# Precipitation
# TODO: check if anything is affected by gaps in precip
flux_qc <- mutate(
  flux_qc,
  QC_LE_PRECIP = if_else(data$P_RAIN > rain_thr, 2, 0),
  QC_FC_PRECIP = if_else(data$P_RAIN > rain_thr, 2, 0)
)

### Low frequency spikes ----

flux_qc$QC_H_PRELIM <- combn_QC(
  as.data.frame(flux_qc), str_subset(names(flux_qc), "_H_"), no_messages = TRUE
)
flux_qc$QC_LE_PRELIM <- combn_QC(
  as.data.frame(flux_qc), str_subset(names(flux_qc), "_LE_"),
  na.as_0_pattern = "PRECIP$", no_messages = TRUE
)
flux_qc$QC_FC_PRELIM <- combn_QC(
  as.data.frame(flux_qc), str_subset(names(flux_qc), "_FC_"),
  na.as_0_pattern = "PRECIP$", no_messages = TRUE
)

flux_qc$timestamp <- data$TIMESTAMP
flux_qc$GR <- data$SW_IN_POT

flux_qc$QC_H_SPIKESLF <- despikeLF(
  as.data.frame(bind_cols(data, flux_qc)), 
  var = "H", qc_flag = "QC_H_PRELIM", 
  var_thr = c(-200, 600), light = "GR", night_thr = night_thr
)

flux_qc$QC_LE_SPIKESLF <- despikeLF(
  as.data.frame(bind_cols(data, flux_qc)), 
  var = "LE", qc_flag = "QC_LE_PRELIM", 
  var_thr = c(-200, 800), light = "GR", night_thr = night_thr
)

flux_qc$QC_FC_SPIKESLF <- despikeLF(
  as.data.frame(bind_cols(data, flux_qc)), 
  var = "FC", qc_flag = "QC_FC_PRELIM", 
  var_thr = c(-50, 50), light = "GR", night_thr = night_thr
)

# FCH4 is checked for range but not spikes
if (site %in% c("JLN", "JLR")) {
  # flux_qc$QC_FCH4_LIM <- if_else(!between(data$FCH4, -0.2, 1.0), 2, 0)
  flux_qc$QC_FCH4_LIM <- if_else(!between(data$FCH4, -0.15, 0.85), 2, 0)
} else if (site == "JLA") {
  # flux_qc$QC_FCH4_LIM <- if_else(!between(data$FCH4, -0.2, 0.5), 2, 0)
  flux_qc$QC_FCH4_LIM <- if_else(!between(data$FCH4, -0.15, 0.3), 2, 0)
}

flux_qc <- select(flux_qc, -timestamp, -GR, -ends_with("_PRELIM"))

### Footprint ----

# source("source/footprint-model.R")
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

## Combine flags and add to dataset --------------------------------------------

met_vars <- unique(str_match(names(met_qc), "QC_(.*)_[^_]*$")[, 2])
met_flags <- met_vars |> 
  set_names(glue("QC_{met_vars}")) |> 
  map(\(x) str_subset(names(met_qc), glue("_{x}_"))) |> 
  map(\(x) combn_QC(as.data.frame(met_qc), x, additive = FALSE, na.as = NA)) |> 
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

# TODO: save these summaries?

## Write quality-flagged dataset with documentation ----------------------------

qc_file <- glue("eddypro_{site}-{year}_fluxnet_adv_qaqc.csv")
write_csv(cleaned, file.path(output_dir, qc_file))

# Documentation
# [goes here]

