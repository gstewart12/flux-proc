# Correct EddyPro output -------------------------------------------------------

# Purpose: 
# Convert units, adjust for instrument offsets, fix known issues

# Load the required packages
library("yaml")
library("httr2")
library("xml2")
library("glue")
library("tidyverse")

## Initialize script settings & documentation ----------------------------------

data_dir <- file.path("data", site) # path to directory containing all data
output_dir <- file.path("output", site, year) # where to write output

# Load site configuration
site_cfg <- read_yaml("config/site-config.yml")[[site]]

# Load metadata file
md <- read_yaml(file.path(data_dir, "metadata.yml"))
lat <- md$tower$coords$lat
lon <- md$tower$coords$lon
tz <- md$tz$utc_offset

## Read merged data file -------------------------------------------------------

data_file <- glue("eddypro_{site}-{year}_fluxnet_adv_merged.csv")

data <- read_csv(file.path(output_dir, data_file))

corr <- data

## Correct variable names & units ----------------------------------------------

# Simplify biomet names
corr <- corr |>
  rename_with(\(x) str_remove(x, "_1_1$")) |>
  rename_with(\(x) str_remove(x, "_1$")) |>
  rename(G_1 = G, SWC_1 = SWC, TS_1 = TS)

# Simplify other names
corr <- corr |>
  # Remove "MEAS" from gas measurement stats
  rename_with(\(x) str_remove(x, "_MEAS")) |>
  # Remove "SINGLE" suffix from storage fluxes (they are all single-level)
  rename_with(\(x) str_remove(x, "_SINGLE")) |>
  # Remove "HF" suffix from random uncertainties
  rename_with(\(x) str_remove(x, "_HF"))

# Specify biomet versions of vars
corr <- rename(corr, TA_MET = TA, RH_MET = RH)

# Calculate version of VPD from biomet vars
corr$VPD_MET <- bigleaf::rH.to.VPD(corr$RH_MET/100, corr$TA_MET) * 10

# Change CH4 units from nmol to umol
ch4_vars <- c(
  "FCH4", "FCH4_RANDUNC", "SCH4", "FCH4_VADV", "CH4_MIXING_RATIO", 
  "CH4", "FCH4_UNCORR", "FCH4_STAGE1", "FCH4_STAGE2"
)
corr <- mutate(corr, across(all_of(ch4_vars), \(x) x/1000))

# Air pressure from Pa to kPa
corr$CUSTOM_AIR_P_MEAN <- corr$CUSTOM_AIR_P_MEAN/1000

# Air temperature from K to C
corr$CUSTOM_AIR_T_MEAN <- corr$CUSTOM_AIR_T_MEAN - 273.15

## Correct fluxes for storage --------------------------------------------------

# Note: this overwrites non-storage corrected fluxes
# TODO: consider filling gaps with non-zero value (e.g. running average)

corr <- mutate(
  corr, 
  H = H + replace_na(SH, 0),
  LE = LE + replace_na(SLE, 0),
  ET = ET + replace_na(SET, 0),
  FH2O = FH2O + replace_na(SH2O, 0),
  FC = FC + replace_na(SC, 0),
  FCH4 = FCH4 + replace_na(SCH4, 0)
)

## Correct wind direction for magnetic declination -----------------------------

# Note: cannot guarantee this will work forever; API may change

# Get declination
noaa_url <- "https://www.ngdc.noaa.gov"
noaa_req <- noaa_url |>
  request() |>
  req_url_path_append(
    "geomag-web","calculators", "calculateUshistoric"
  ) |>
  req_url_query(
    lat1 = lat, lon1 = lon, 
    startYear = year, endYear = year, 
    key = "yKc9F", resultFormat = "xml"
  ) 
mag_decl <- noaa_req |>
  req_perform() |> 
  resp_body_xml() |> 
  xml_find_first(".//declination") |> 
  xml_double()

# Reassign magnetic wind direction; apply geographic offset
corr <- mutate(
  corr,
  WD_MAG = WD,
  # See https://www.nwcg.gov/course/ffm/location/65-declination
  # - explains why declination is subtracted (i.e. negative value is added)
  WD = (WD + mag_decl) %% 360,
  .after = WD
)

## Correct incoming radiation --------------------------------------------------

# Re-calculate potential radiation
corr <- corr |>
  mutate(
    SW_IN_POT = bigleaf::potential.radiation(
      yday(TIMESTAMP), hour(TIMESTAMP) + minute(TIMESTAMP)/60, lat, lon, tz
    )
  )

# Correct for zero-offset
offset <- corr |>
  filter(SW_IN_POT == 0) |>
  summarize(across(c(SW_IN, PPFD_IN), \(x) mean(x, na.rm = TRUE))) |>
  as.list()
corr <- corr |>
  # Adjust all data (day & night) by offset value
  mutate(
    SW_IN = SW_IN - offset$SW_IN,
    PPFD_IN = PPFD_IN - offset$PPFD_IN
  ) |>
  # Set all nighttime values to 0
  # mutate(across(c(SW_IN, PPFD_IN), \(x) if_else(SW_IN_POT == 0, 0, x))) |>
  mutate(
    # Set negative values to 0 if potential radiation is 0
    across(c(SW_IN, PPFD_IN), \(x) if_else(x < 0 & SW_IN_POT == 0, 0, x)),
    # Set remaining negative values to NA
    across(c(SW_IN, PPFD_IN), \(x) if_else(x < 0, NA, x))
  )

## Correct known site/year-specific issues -------------------------------------

# Get excluded periods from file
# Note: DO NOT edit these files in Excel - it messes with date formatting
# TODO: make code flexible enough to allow date format to change
# exclude <- as.list(read_csv(file.path(data_dir, "flagged_periods.csv")))
exclude <- site_cfg$flagged_periods
exclude_log <- c()

# Loop through vars and remove data within specified time ranges
# for (i in 1:length(exclude$var)) {
#   corr <- mutate(
#     corr, 
#     "{exclude$var[i]}" := if_else(
#       between(TIMESTAMP, exclude$from[i], exclude$to[i]), 
#       NA_real_, .data[[exclude$var[i]]]
#     )
#   )
# }
for (var in names(exclude)) {
  for (i in 1:length(exclude[[var]])) {
    period <- exclude[[var]][[i]]
    
    # Get flagged bounds within processing year
    from <- max(ymd_hm(period$from), first(corr$TIMESTAMP))
    to <- min(ymd_hm(period$to), last(corr$TIMESTAMP))
    
    # Skip if nothing flagged within processing year
    if (from > to) next
    
    # Remove data between flagged bounds (inclusive)
    if (!is.null(period$remove)) {
      # Remove selected data
      var_expr <- paste(var, period$remove)
      var_flag <- eval(parse(text = var_expr), envir = corr)
      corr[[var]] <- if_else(
        between(corr$TIMESTAMP, from, to) & var_flag,
        NA_real_, corr[[var]]
      )
    } else {
      # Remove all data
      var_expr <- var
      corr[[var]] <- if_else(
        between(corr$TIMESTAMP, from, to),
        NA_real_, corr[[var]]
      )
    }
    
    # Write to log for documentation
    exclude_log <- c(
      exclude_log,
      glue("  - {var_expr} from {from} to {to}; reason: {period$reason}")
    )
  }
}

# Propagate missing data in biomet logger power
# - 'LOGGERPOWER' is used as a proxy to flag all biomet data
corr <- mutate(
  corr, across(LOGGERTEMP:TS_3, \(x) if_else(is.na(LOGGERPOWER), NA, x))
)

## Set primary versions of replicate measurements ------------------------------

# Set PA_EP to NA when missing CH4 data (PA_EP is from LI-7700)
corr$PA_EP <- if_else(
  is.na(corr$CH4) & !is.na(data$CUSTOM_AIR_P_MEAN), NA, corr$PA_EP
)

# Since LI-7700 TA & PA are flaky, is a better solution to run EddyPro with 
# those vars from LI-7500? --> tested this: no, not really

# Get site-specific configuration
# rep_vars <- config$var_merge
# 
# for (var in names(rep_vars)) {
#   method <- rep_vars[[var]]$method
#   reps <- corr[rep_vars[[var]]$order]
#   
#   if (method == "coalesce") {
#     corr[[var]] <- coalesce(!!!as.list(reps))
#   } else if (method == "average") {
#     corr[[var]] <- rowMeans(cbind(reps), na.rm = TRUE)
#   }
# }

# TODO: put everything here in metadata files, use above commented code

if (site == "JLA") {
  
  corr$TA <- coalesce(corr$CUSTOM_AIR_T_MEAN, corr$TA_EP, corr$TA_MET)
  corr$PA <- coalesce(corr$CUSTOM_AIR_P_MEAN, corr$PA_EP)
  corr$VPD <- coalesce(corr$VPD_EP, corr$VPD_MET)
  corr$RH <- coalesce(corr$RH_EP, corr$RH_MET)
  
  # Average reps
  corr$G <- rowMeans(cbind(corr$G_1, corr$G_2, corr$G_3))
  corr$TS <- rowMeans(cbind(corr$TS_1, corr$TS_2, corr$TS_3))
  corr$SWC <- rowMeans(cbind(corr$SWC_1, corr$SWC_2, corr$SWC_3))
  
  # Remove bias due to varying rep availability
  # Approximation of all reps from just reps 1 & 3
  corr$TS <- coalesce(
    corr$TS, rowMeans(cbind(corr$TS_1, corr$TS_3)) * 0.9973 + 0.0716
  )
  corr$SWC <- coalesce(
    corr$SWC, rowMeans(cbind(corr$SWC_1, corr$SWC_3)) * 0.8918 + 3.4879
  )
  
} else if (site == "JLN") {
  
  corr$TA <- coalesce(corr$CUSTOM_AIR_T_MEAN, corr$TA_EP, corr$TA_MET)
  corr$PA <- coalesce(corr$CUSTOM_AIR_P_MEAN, corr$PA_EP)
  corr$VPD <- coalesce(corr$VPD_EP, corr$VPD_MET)
  corr$RH <- coalesce(corr$RH_EP, corr$RH_MET)
  
  # Average reps
  # G_2 inverted beginning in Jun 2019
  corr$G_2 <- if_else(
    corr$TIMESTAMP >= ymd_hm("2019-06-07 11:45"), -corr$G_2, corr$G_2
  )
  corr$G <- rowMeans(cbind(corr$G_1, corr$G_2, corr$G_3))
  corr$TS <- rowMeans(cbind(corr$TS_1, corr$TS_2, corr$TS_3))
  corr$SWC <- rowMeans(cbind(corr$SWC_1, corr$SWC_2, corr$SWC_3))
  
  # Remove bias due to varying rep availability
  # Approximation of all reps from just reps 1 & 3
  corr$TS <- coalesce(
    corr$TS, rowMeans(cbind(corr$TS_1, corr$TS_3)) * 1.0111 + -0.1904
  )
  corr$SWC <- coalesce(
    corr$SWC, rowMeans(cbind(corr$SWC_1, corr$SWC_3)) * 0.8061 + 12.5992
  )
  corr$G <- coalesce(
    corr$G, rowMeans(cbind(corr$G_1, corr$G_3)) * 0.43043 + -0.03457
  )
  # Approximation of all reps from just rep 3
  corr$TS <- coalesce(corr$TS, corr$TS_3 * 1.0189 + -0.2274)
  corr$SWC <- coalesce(corr$SWC, corr$SWC_3 * 0.8359 + 10.2531)
  
} else if (site == "JLR") {
  
  corr$TA <- coalesce(corr$TA_MET, corr$CUSTOM_AIR_T_MEAN, corr$TA_EP)
  corr$PA <- coalesce(corr$CUSTOM_AIR_P_MEAN, corr$PA_EP)
  corr$VPD <- coalesce(corr$VPD_MET, corr$VPD_EP)
  corr$RH <- coalesce(corr$RH_MET, corr$RH_EP)
  
  # Average reps
  # G_1 inverted throughout
  corr$G_1 <- -corr$G_1
  corr$G <- rowMeans(cbind(corr$G_1, corr$G_2, corr$G_3))
  # TS_1 has WAY more variance than other reps - maybe it came loose...
  corr$TS <- rowMeans(cbind(corr$TS_1, corr$TS_2, corr$TS_3))
  corr$SWC <- rowMeans(cbind(corr$SWC_1, corr$SWC_2, corr$SWC_3))
  
  # Remove bias due to varying rep availability
  # Approximation of all reps from just reps 1 & 3
  corr$TS <- coalesce(
    corr$TS, rowMeans(cbind(corr$TS_1, corr$TS_3)) * 0.9882 + 0.2767
  )
  corr$SWC <- coalesce(
    corr$SWC, rowMeans(cbind(corr$SWC_1, corr$SWC_3)) * 0.8124 + 11.7965
  )
  # Approximation of all reps from just rep 1
  corr$TS <- coalesce(corr$TS, corr$TS_1 * 0.825 + 2.670)
  corr$SWC <- coalesce(corr$SWC, corr$SWC_1 * 0.4142 + 37.9382)
  
}

# Calculate water vapor partial pressure at saturation
calc_esat <- function(x) {
  # x is air temperature in degC
  x <- x + 273.15 # degC to degK
  esat_pa <- x^(-8.2) * exp(77.345 + 0.0057 * x - 7235/x)
  esat_pa/100 # units = hPa
}

# Calculate dew point temperature
calc_tdew <- function(x) {
  # x is water vapor partial pressure in hPa
  x <- x/10 # hPa to kPa
  240.97 * (log(x/0.611)/(17.502 - log(x/0.611))) # units = degC
}

# Update atmospheric conditions
corr$TDEW_EP <- corr$TDEW
corr$VAPOR_PARTIAL_PRESSURE_SAT <- calc_esat(corr$TA)
corr$VAPOR_PARTIAL_PRESSURE <- corr$RH/100 * corr$VAPOR_PARTIAL_PRESSURE_SAT
# corr$TDEW <- bigleaf::dew.point(corr$TA, corr$VPD * 10) # too slow
corr$TDEW <- calc_tdew(corr$VAPOR_PARTIAL_PRESSURE)
# TODO also update VPD, RH?


## Write corrected dataset with documentation ----------------------------------

corr_file <- glue("eddypro_{site}-{year}_fluxnet_adv_corrected.csv")
write_csv(corr, file.path(output_dir, corr_file))

# Documentation
write_lines(
  c(
    "EddyPro output files corrected via '02-correct.R'",
    "",
    "Corrections",
    glue(
      "  - Converted units from nmol to umol in vars: ", 
           glue_collapse(ch4_vars, ", ")
    ),
    glue("  - Converted 'CUSTOM_AIR_P_MEAN' from Pa to kPa"),
    glue("  - Converted 'CUSTOM_AIR_T_MEAN' from K to C"),
    glue("  - Applied storage corrections to all fluxes"),
    glue("  - Assigned var 'WD' to new var 'WD_MAG'"),
    glue("  - Corrected 'WD' for magnetic declination ({mag_decl} in {year})"),
    glue("  - Recalculated 'SW_IN_POT' using bigleaf::potential.radiation"),
    glue(
      "  - Corrected zero-offsets in 'SW_IN' ({offset$SW_IN})", 
         " and 'PPFD_IN' ({offset$PPFD_IN})"
    ),
    glue("  - Set negative nighttime values of 'SW_IN' and 'PPFD_IN' to zero"),
    "",
    "Data removed",
    glue_collapse(exclude_log, "\n"),
    "",
    glue("File out: {corr_file}"),
    "",
    "Session info",
    glue(str_pad("date", 10, "right"), as.character(today())),
    glue(str_pad("version", 10, "right"), R.version.string),
    glue(str_pad("os", 10, "right"), osVersion)
  ),
  file.path(output_dir, str_replace(corr_file, "\\.csv", "\\.txt"))
)
