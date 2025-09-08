# Gap-fill non-CH4 flux variables ----------------------------------------------

# Purpose: Fill all gaps in non-CH4 fluxes using marginal distribution sampling, 
# partition NEE into Reco & GPP using both nighttime & daytime algorithms.

# Load the required packages
library("openeddy")
library("REddyProc")
library("yaml")
library("glue")
library("tidyverse")

## Initialize script settings & documentation ----------------------------------

site <- "JLR" # three-letter site code
years <- 2019:2023 # years to process
output_dirs <- file.path("output", site, years) # where to write output
# Met-gapfilled data
data_files <- glue("eddypro_{site}-{years}_fluxnet_adv_metgf.csv")
datas <- read_csv(file.path(output_dirs, data_files))

# Set nighttime threshold for incoming radiation
# - REddyProc default: night = Rg_orig <= 10 & PotRad_NEW == 0
# - lasslop2010 use night = Rg_orig <= 20
night_thr <- 10

# Set u* filtering options
ustar_thr <- 0.05 # nominal threshold
flag_next <- FALSE # whether to also flag next record when ustar < ustar_thr
filter_daytime <- TRUE # whether to apply ustar_thr to day as well as night

## Perform yearly MDS gap-filling

for (year in years) {
  
  output_dir <- file.path("output", site, year) # where to write output
  
  # Load metadata file
  md <- read_yaml(file.path("data", site, "metadata.yml"))
  
  ## Read required data files --------------------------------------------------
  
  # Met-gapfilled data
  data_file <- glue("eddypro_{site}-{year}_fluxnet_adv_metgf.csv")
  
  data <- read_csv(file.path(output_dir, data_file))
  
  # Apply QC flags
  cleaned <- mutate(
    data,
    H = apply_QC(H, QC_H),
    LE = apply_QC(LE, QC_LE),
    FC = apply_QC(FC, QC_FC)
  )
  
  # Format data for REddyProc
  data_mds <- cleaned |>
    mutate(TIMESTAMP = TIMESTAMP + minutes(15)) |>
    select(
      timestamp = TIMESTAMP,
      H, LE, NEE = FC, Rg = SW_IN_F, Tair = TA_F, VPD = VPD_F, Ustar = USTAR_F
    ) |>
    as.data.frame()
  
  # Gap-filling with & without u* threshold ------------------------------------
  
  # Initialize REddyProc class
  proc <- sEddyProc$new(
    ID = site, 
    Data = data_mds, 
    ColNames = c("H", "LE", "NEE", "Rg", "Tair", "VPD", "Ustar"), 
    ColPOSIXTime = "timestamp"
  )
  
  # Set location info
  proc$sSetLocationInfo(
    LatDeg = md$tower$coords$lat, 
    LongDeg = md$tower$coords$lon, 
    TimeZoneHour = md$tz$utc_offset
  )
  
  # Fill gaps in energy fluxes with MDS
  proc$sMDSGapFill("H")
  proc$sMDSGapFill("LE")
  
  # Estimate u* threshold
  # proc$sEstUstarThold()

  # Fill gaps in NEE with MDS
  proc$sMDSGapFill("NEE")
  # Separate versions for with/without u* threshold
  proc$sMDSGapFillAfterUstar(
    "NEE", 
    uStarTh = ustar_thr, 
    isFilterDayTime = filter_daytime,
    isFlagEntryAfterLowTurbulence = flag_next
  )
  
  # Fill gaps in met vars (shouldn't actually be any gaps)
  proc$sMDSGapFill("Rg", FillAll = FALSE)
  proc$sMDSGapFill("Tair", FillAll = FALSE)
  proc$sMDSGapFill("VPD", FillAll = FALSE)
  
  # Partition NEE (nighttime-based algorithm)
  proc$sMRFluxPartition("NEE_f")
  proc$sMRFluxPartition("NEE_f", suffix = "uStar")
  
  # Partition NEE (daytime-based algorithm)
  proc$sGLFluxPartition("NEE_f")
  proc$sGLFluxPartition("NEE_f", suffix = "uStar")
  
  # Note: in REddyProc, night = Rg_orig <= 10 & PotRad_NEW == 0
  
  # Gather results
  results <- proc$sExportResults()
  mds_filled <- select(
    results, 
    H_F = H_f, H_F_QC = H_fqc, LE_F = LE_f, LE_F_QC = LE_fqc, 
    NEE_F = NEE_f, NEE_F_QC = NEE_fqc, 
    NEE_USTAR_F = NEE_uStar_f, NEE_USTAR_F_QC = NEE_uStar_fqc
  )
  nee_part <- select(
    results, 
    RECO_NT = Reco, GPP_NT = GPP_f, RECO_DT = Reco_DT, GPP_DT, 
    RECO_NT_USTAR = Reco_uStar, GPP_NT_USTAR = GPP_uStar_f, 
    RECO_DT_USTAR = Reco_DT_uStar, GPP_DT_USTAR = GPP_DT_uStar
  )
  
  processed <- bind_cols(data, mds_filled, nee_part)
  
  # Quick view of gap-filled variables
  processed |>
    ggplot(aes(TIMESTAMP, H_F)) +
    geom_point(size = 0.75)
  processed |>
    ggplot(aes(TIMESTAMP, LE_F)) +
    geom_point(size = 0.75)
  processed |>
    mutate(QC = factor(NEE_F_QC)) |>
    ggplot(aes(TIMESTAMP, NEE_F, color = QC)) +
    geom_point(size = 0.75)
  processed |>
    mutate(QC = factor(NEE_USTAR_F_QC)) |>
    ggplot(aes(TIMESTAMP, NEE_USTAR_F, color = QC)) +
    geom_point(size = 0.75)
  processed |>
    ggplot(aes(TIMESTAMP, GPP_NT)) +
    geom_point(size = 0.75)
  processed |>
    ggplot(aes(TIMESTAMP, GPP_DT)) +
    geom_point(size = 0.75)
  processed |>
    ggplot(aes(TIMESTAMP, GPP_NT_USTAR)) +
    geom_point(size = 0.75)
  processed |>
    ggplot(aes(TIMESTAMP, GPP_DT_USTAR)) +
    geom_point(size = 0.75)
  
  # Annual CO2 balance (g C m-2)
  # summarize(
  #   processed, 
  #   across(c(NEE_F, NEE_USTAR_F, RECO_NT:GPP_DT_USTAR), \(x) sum(x) * 0.021618)
  # )
  
  # Estimate annual budget & uncertainty ---------------------------------------
  
  budget <- results |>
    mutate(timestamp = processed$TIMESTAMP) |>
    agg_sum("%Y", agg_per = "y-1", NEE_scor = FALSE) |>
    select(
      Intervals, days, contains("_f_"), Reco_sum, Reco_uStar_sum,
      Reco_DT_sum, GPP_DT_sum, Reco_DT_uStar_sum, GPP_DT_uStar_sum,
      -starts_with("Rg"), -starts_with("Tair"), -starts_with("VPD")
    ) |>
    mutate(
      NEP_DT_sum = Reco_DT_sum - GPP_DT_sum,
      NEP_DT_uStar_sum = Reco_DT_uStar_sum - GPP_DT_uStar_sum
    )
  
  # Aggregate uncertainty
  budget_unc <- results |>
    mutate(timestamp = processed$TIMESTAMP) |>
    agg_fsd("%Y", agg_per = "y-1") |>
    pluck("sum") |>
    select(where(\(x) !all(is.na(x))))
  
  # Reco and GPP uncertainty
  budget_unc_dt <- results |>
    mutate(timestamp = processed$TIMESTAMP) |>
    select(-Reco_DT, -Reco_DT_SD, -GPP_DT, -GPP_DT_SD) |>
    agg_DT_SD("%Y", agg_per = "y-1") |>
    pluck("sum")
  
  # Note:
  # Uncertainty products of agg_fsd and agg_DT_SD are reported as std. devs.
  # To represent uncertainty bounds for given confidence interval:
  # e.g. SD * 1.96 for 95% confidence level
  
  # Budgets estimated with nighttime & daytime methods
  aggregated <- budget |>
    left_join(budget_unc) |>
    left_join(budget_unc_dt)
  
  agg_file <- glue("mds-budgets-{site}-{year}.csv")
  write_csv(aggregated, file.path(output_dir, agg_file))
  
  ## Write gap-filled dataset with documentation -------------------------------
  
  proc_file <- glue("eddypro_{site}-{year}_fluxnet_adv_processed.csv")
  write_csv(processed, file.path(output_dir, proc_file))
  
  # Documentation
  # [goes here]
  
}

# Summarize all yearly budgets (g C m-2)
file.path(
  "output", site, years, 
  glue("eddypro_{site}-{years}_fluxnet_adv_processed.csv")
) |>
  read_csv() |>
  group_by(year = year(TIMESTAMP)) |>
  summarize(
    across(c(NEE_F, NEE_USTAR_F, RECO_NT:GPP_DT_USTAR), \(x) sum(x) * 0.021618)
  ) |>
  mutate(
    NEE_DT = RECO_DT - GPP_DT,
    NEE_DT_USTAR = RECO_DT_USTAR - GPP_DT_USTAR,
    .after = NEE_USTAR_F
  )


# Gap-filling with u* threshold & uncertainty ----------------------------------

# plot(NEE ~ timestamp, data_mds)
# # Estimate distribution of u* thresholds
# proc$sEstUstarThold(
#   ctrlUstarEst = usControlUstarEst(),
#   ctrlUstarSub = usControlUstarSubsetting(),
# )
# # proc$sEstUstarThreshold()
# proc$sEstimateUstarScenarios()
# # Check results
# proc$sGetEstimatedUstarThresholdDistribution()
# proc$sGetUstarScenarios()
# proc$sPlotNEEVersusUStarForSeason()

# # Fill gaps in NEE with MDS
# proc$sMDSGapFillUStarScens("NEE")
# 
# # Fill gaps in met vars
# proc$sMDSGapFill("Tair", FillAll = FALSE, minNWarnRunLength = NA)
# proc$sMDSGapFill("VPD", FillAll = FALSE, minNWarnRunLength = NA)
# 
# # Partition NEE (nighttime-based algorithm)
# proc$sMRFluxPartitionUStarScens()
# 
# # Partition NEE (daytime-based algorithm)
# proc$sGLFluxPartitionUStarScens()

# REddyProc variable & suffix definitions --------------------------------------

# Gap-filling output
# *_orig - 
# *_f - 
# *_fqc - 
# *_fall - 
# *_fall_qc - 
# *_fnum - 
# *_fsd - 
# *_fmeth - 
# *_fwin - 

# Nighttime-based partitioning output
# PotRad - Potential radiation
# FP_NEEnight - Good (original) NEE nighttime fluxes used for flux partitioning
# FP_Temp - Good (original) temperature measurements used for flux partitioning
# E_0 - Estimated temperature sensitivity
# R_ref - Estimated reference respiration
# Reco - Estimated ecosystem respiration
# GPP_f - Estimated gross primary production

# Daytime-based partitioning output
# Reco_DT -	predicted ecosystem respiration
# GPP_DT - predicted gross primary production
# FP_VARnight - NEE filtered for nighttime records (others NA)
# FP_VARday - NEE filtered for daytime records (others NA)
# NEW_FP_Temp - temperature after filtering for quality flag
# NEW_FP_VPD - VPD after filtering for quality flag
# FP_RRef_Night - basal respiration estimated from nighttime
# FP_qc - quality flag: 
#   0: good parameter fit 
#   1: some parameters out of range, required refit
#   2: next parameter estimate is more than two weeks away
# FP_dRecPar - records until/after closest record with a parameter estimate
# FP_errorcode - why LRC-fit was not successful or was rejected
# FP_GPP2000 - predicted GPP at VPD = 0 and PAR = 2000
#   - a surrogate for maximum photosynthetic capacity


