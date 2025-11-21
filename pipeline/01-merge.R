# Merge EddyPro output by year -------------------------------------------------

# Purpose: 
# Combine output from different EddyPro runs, check formatting, add external 
# data, save a single file ready for post-processing

# Load the required packages
library("yaml")
library("glue")
library("readxl")
library("tidyverse")

# Load custom functions
source("R/pipeline-funs.R")

## Initialize script settings & documentation ----------------------------------

data_dir <- file.path("data", site) # path to directory containing all data
output_dir <- file.path("output", site, year) # where to write output

# If first time running site/year, need to create output directory
if (!dir.exists(output_dir)) dir.create(output_dir)

# Load configuration file
cfg <- read_yaml("config/config.yml")

## Read and combine the EddyPro files ------------------------------------------

data_files <- get_eddypro_files_year(data_dir, year)
biomet_files <- get_eddypro_files_year(data_dir, year, type = "biomet")

# Read data files
skip_cols <- c("VERSION_1_1_1", "_0_0_1") # these cause errors and have no data
na_vals <- cfg$merge$na_vals # read these as NA
data <- data_files |>
  map(\(x) read_csv(x, na = na_vals, col_select = -any_of(skip_cols))) |>
  bind_rows() |>
  # Subset data from selected year
  filter(str_sub(TIMESTAMP_START, 1, 4) == year) |>
  arrange(TIMESTAMP_START)
# Note: parsing problems come from skipped cols and are safe to ignore

# Handle edge case - two versions of the same column (duplicate names)
dupe_cols <- str_subset(names(data), "\\.{3}")
if (length(dupe_cols)) {
  dupe_col_stems <- str_remove(dupe_cols, "\\.{3}.+")
  for (col in unique(dupe_col_stems)) {
    dupe_col_data <- data |>
      select(all_of(c(col, dupe_cols[str_starts(dupe_cols, col)]))) |>
      as.list()
    # Merge duplicates with original
    data[[col]] <- coalesce(!!!dupe_col_data)
    # Remove duplicate versions
    data <- select(data, -all_of(dupe_cols[str_starts(dupe_cols, col)]))
  }
}

merged <- data

## Patch biomet data -----------------------------------------------------------

### Clear erroneous data ----

# Undo LOCF (last observation carried forward) applied by EddyPro
# - seems to only affect fluxnet output file

# Read biomet output file (no LOCF)
biomet <- biomet_files |>
  map(\(x) read_csv(x, na = na_vals)) |>
  map(\(x) slice(x, -1)) |>
  bind_rows() |>
  mutate(TIMESTAMP = ymd_hm(paste(date, time)) - minutes(15), .before = 1) |>
  filter(year(TIMESTAMP) == year) |>
  arrange(TIMESTAMP) |>
  type_convert()

merged <- merged |>
  mutate(
    # Overwrite biomet logger power
    LOGGERPOWER_1_1_1 = as.numeric(biomet$LOGGERPOWER_1_1_1),
    # Propagate NAs to rest of biomet vars
    across(
      c(LOGGERTEMP_1_1_1:last_col()), 
      \(x) if_else(is.na(LOGGERPOWER_1_1_1), NA_real_, x)
    ),
    # Remove nighttime determinations based on LOCF biomet data
    NIGHT = if_else(is.na(SW_IN_1_1_1), NA_real_, NIGHT)
  )

### Add any recovered data ----

# TODO: need to rethink this whole step
# - for biomet, could just process it all separately

# Check if any recovered summary files exist
rec_data_files <- list.files(
  dirname(dirname(data_files)), pattern = "data_recovered", full.names = TRUE
)
# TEMP fix for recovered data but no eddypro output
# TODO: real fix
if (identical(site, "JLN") & year == 2023) {
  rec_data_files <- c(
    rec_data_files, 
    file.path(data_dir, "1773_2023-10-02/data_recovered_1773_2023-10-02.csv")
  )
}
if (length(rec_data_files) > 0) {
  
  rec_data <- rec_data_files |>
    read_csv() |>
    select(-contains("SHFSENS")) |>
    rename_with(\(x) str_replace(x, "SHF", "G")) |>
    rename_with(\(x) str_replace(x, "LW", "LW_")) |>
    rename_with(\(x) str_replace(x, "SW", "SW_"), -starts_with("SWC")) |>
    rename_with(\(x) str_replace(x, "PPFD_1", "PPFD_IN_1")) |>
    rename_with(\(x) str_replace(x, "RN_1", "NETRAD_1")) |>
    # Summaries have EddyPro units, so need to convert to Fluxnet
    mutate(
      across(starts_with("SWC"), \(x) x * 100),
      across(c(T_SONIC, TA_EP, TDEW), \(x) x - 273.15),
      across(c(starts_with("TA_1"), starts_with("TS_")), \(x) x - 273.15),
      PA_EP = PA_EP/1000,
      across(c(VPD_EP, starts_with("VAPOR_PARTIAL")), \(x) x/100),
      across(starts_with("P_RAIN"), \(x) x * 1000)
    )
  # Recovered because no GHG files available, so these are new rows
  merged <- rows_insert(
    merged, rec_data, by = c("TIMESTAMP_START", "TIMESTAMP_END"), 
    conflict = "ignore"
  )
}

# Check if any recovered biomet files exist
rec_biomet_files <- list.files(
  dirname(dirname(data_files)), pattern = "biomet_recovered", full.names = TRUE
)
# TEMP fix for recovered data but no eddypro output
# TODO: real fix
if (identical(site, "JLR") & year == 2023) {
  rec_biomet_files <- c(
    rec_biomet_files, 
    file.path(data_dir, "1772_2023-10-03/biomet_recovered_1772_2023-10-03.csv")
  )
}
if (length(rec_biomet_files) > 0) {
  
  rec_biomet <- rec_biomet_files |>
    read_csv() |>
    select(-contains("SHFSENS")) |>
    rename_with(\(x) str_replace(x, "SHF", "G")) |>
    rename_with(\(x) str_replace(x, "LW", "LW_")) |>
    rename_with(\(x) str_replace(x, "SW", "SW_"), -starts_with("SWC")) |>
    rename(PPFD_IN_1_1_1 = PPFD_1_1_1, NETRAD_1_1_1 = RN_1_1_1) |>
    mutate(across(starts_with("SWC"), \(x) x * 100))
  
  # Recovered because EddyPro skipped GHG files, so these rows are patched
  merged <- rows_patch(
    merged, rec_biomet, by = c("TIMESTAMP_START", "TIMESTAMP_END"), 
    unmatched = "ignore"
  )
  
}

## Apply formatting corrections to combined dataset ----------------------------

# Remove empty/redundant/irrelevant columns 
merged <- select(
  merged,
  -contains("NONE"), -contains("LI7200"), -contains("TUBE"), -contains("CELL"), 
  -contains("H_BU"), -contains("KRYPTON"), -contains("IBROM"),
  -contains("CH4_TC"), -contains("SHFSENS"), -contains("CUSTOM_CH4_AUX"),
  -contains("CUSTOM_AUXILIARY_INPUT"),
  -any_of(
    c(
      "INST_LI7700_RSSI", "ROT_ROLL", 
      "CUSTOM_U_MEAN", "CUSTOM_V_MEAN", "CUSTOM_W_MEAN", "CUSTOM_TS_MEAN"
    )
  )
)
dropped_cols <- c(setdiff(names(data), names(merged)), tolower(skip_cols))

# Make timestamp, ensure it forms a regular sequence for the entire year
year_seq <- seq(
  ymd_hm(glue("{year}-01-01 00:15")), ymd_hm(glue("{year}-12-31 23:45")), 
  by = "30 mins"
)
merged <- merged |>
  # Center timestamp in middle of averaging period
  mutate(TIMESTAMP = ymd_hm(TIMESTAMP_START) + minutes(15), .before = 1) |>
  right_join(as_tibble_col(year_seq, "TIMESTAMP"), by = join_by(TIMESTAMP)) %>%
  # Select duplicate rows with least amount of missing data
  arrange(TIMESTAMP, rowSums(is.na(.))) |> 
  distinct(TIMESTAMP, .keep_all = TRUE)
n_empty <- nrow(merged) - nrow(data)

## Join external data sources --------------------------------------------------

### ERA variables ----

era5_file <- glue("output/ERA5/era5-sl-{year}-hh.csv")

# Boundary layer height: necessary for footprint modeling
blh <- era5_file |>
  read_csv() |>
  select(TIMESTAMP, BLH_ERA = BLH)

merged <- left_join(merged, blh)

## Write merged files with documentation ---------------------------------------

# Full dataset
merged_file <- glue("eddypro_{site}-{year}_fluxnet_adv_merged.csv")
write_csv(merged, file.path(output_dir, merged_file))

# TODO: Subset essential variables (required for processing/analysis)
# essentials <- select(
#   merged, ...
# )
# essentials_file <- glue("eddypro_{site}-{year}_fluxnet_adv_merged_essentials.csv")
# write_csv(essentials, file.path(output_dir, essentials_file))

# Documentation
write_lines(
  c(
    "EddyPro output files merged via '01-merge.R'",
    "",
    "EddyPro file(s) merged",
    glue("  - {data_files}"),
    "",
    "Corrections",
    glue("  - Values set to NA: ", glue_collapse(na_vals, ", ")),
    glue("  - Added empty rows (n={n_empty}) to complete full-year sequence"),
    glue("  - Removed columns: ", glue_collapse(dropped_cols, ", ")),
    "",
    glue(
      "Added ERA5 variables: ", 
      glue_collapse(str_subset(names(merged), "_era$"), ", ")
    ),
    "",
    glue("File out: {merged_file}"),
    "",
    "Session info",
    glue(str_pad("date", 10, "right"), as.character(today())),
    glue(str_pad("version", 10, "right"), R.version.string),
    glue(str_pad("os", 10, "right"), osVersion)
  ),
  file.path(output_dir, str_replace(merged_file, "\\.csv", "\\.txt"))
)
