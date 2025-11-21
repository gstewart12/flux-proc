# Run pipeline for selected site(s) and year(s)

# This runs steps 1-3 by site/year, then steps 4-5 by site

# sites <- c("JLA", "JLN", "JLR")
# years <- 2019:2025
# steps <- 1:3

sites <- readline(
  prompt = "Enter site code(s) (JLA, JLN, and/or JLR) separated by commas, or press Enter for all sites: "
)

if (!nchar(sites)) {
  sites <- c("JLA", "JLN", "JLR")
} else {
  # Split, strip, and make uppercase
  sites <- strsplit(sites, ",")[[1]]
  sites <- sites |> trimws() |> toupper() |> unique()
  if (any(!sites %in% c("JLA", "JLN", "JLR"))) {
    stop("Invalid site(s) entered, must be JLA, JLN, and/or JLR.")
  }
}


years <- readline(
  prompt = "Enter year(s) separated by commas, or press Enter for all years: "
)

if (!nchar(years)) {
  years <- 2019:2025
} else {
  # Split, strip, and convert to integer
  years <- strsplit(years, ",")[[1]]
  years <- years |> trimws() |> unique()
  if (any(!years %in% 2019:2025)) {
    stop("Invalid year(s) entered, must be between 2019 and 2025.")
  }
  years <- as.integer(years)
}

steps <- readline(
  prompt = "Enter step(s) separated by commas, or press Enter for all steps: "
)

if (!nchar(steps)) {
  steps <- 1:5
} else {
  # Split, strip, and convert to integer
  steps <- strsplit(steps, ",")[[1]]
  steps <- steps |> trimws() |> unique()
  if (any(!steps %in% 1:5)) {
    stop("Invalid step(s) entered, must be number(s) 1 through 5.")
  }
  steps <- as.integer(steps)
}

for (site in sites) {
  for (year in years) {
    if (1 %in% steps) {
      source("pipeline/01-merge.R")
    }
    if (2 %in% steps) {
      source("pipeline/02-correct.R")
    }
    if (3 %in% steps) {
      source("pipeline/03-qc-auto.R")
    }
  }
}

for (site in sites) {
  if (4 %in% steps) {
    source("pipeline/04-biomet-gapfill.R")
  }
  if (5 %in% steps) {
    source("pipeline/05-mds-gapfill.R")
  }
}
