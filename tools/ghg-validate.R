
# This script checks raw data files (.ghg) for potential issues that can cause 
# errors while running EddyPro

library("tidyverse")

# Functions that should be loaded from a source file

ghg_has_biomet <- function(file) {
  
  contents <- purrr::safely(unzip)(file, list = TRUE)
  
  if (!is.null(contents$error)) {
    return(FALSE)
  } else {
    name <- contents$result$Name
  }
  
  any(stringr::str_detect(name, "biomet.data"))
}

ghg_check_biomet <- function(file) {
  
  nl <- 5
  
  # Check if there is a biomet file
  biomet_exists <- ghg_has_biomet(file)
  data_spec <- readr::cols_only()
  
  if (biomet_exists) {
    name <- basename(file)
    
    data_spec <- readr::spec_tsv(
      unz(file, stringr::str_replace(name, ".ghg", "-biomet.data")), 
      col_types = readr::cols(), 
      na = c("", "NA", "-9999"), 
      skip = nl, 
      progress = FALSE
    )
  }
  
  data_spec
}

# Returns the averaging period (30 min) that a ghg file will be included in
get_time_from_filename <- function(file, round = FALSE) {
  
  name <- basename(file)
  time <- lubridate::ymd_hms(stringr::str_sub(name, 1, 17))
  
  if (round) {
    time <- lubridate::ceiling_date(time, "30 mins", change_on_boundary = TRUE)
  }
  
  time
}

# Set path to raw data directory
dir <- "/Volumes/FLUXDATA/1772_JLR/1772_2024-01-12/raw"
files <- list.files(dir, full.names = TRUE, recursive = TRUE)

table(fs::path_ext(files))

# Biomet file validation
headers <- map(files, safely(ghg_check_biomet), .progress = TRUE)
check <- tibble(
  file = files,
  avg_per = get_time_from_filename(files, round = TRUE),
  ncol = headers |>
    map("result") |>
    map("cols") |>
    map_dbl(length)
)

# Check for ghg files that cannot be extracted
check |> 
  filter(map_lgl(headers, \(x) is.null(x$result)))

# Check for unusual timestamps
check |>
  filter(avg_per > now() | year(avg_per) < 2017)

# Check for unusual number of columns
# - **these files must be archived**
check |>
  filter(ncol < 29, ncol > 0)

# Check for column number change during averaging period
# - **these files must be archived**
check |>
  filter(ncol > 0, ncol != 5) |>
  group_by(avg_per) |>
  filter(n() > 1) |>
  mutate(ncol_diff = ncol - first(ncol)) |>
  filter(ncol_diff != 0)

# Move files to archive
# check |>
#   filter(ncol < 29, ncol > 0) |>
#   pull(file) |>
#   fs::file_move(str_replace(dir, "raw", "archive"))

