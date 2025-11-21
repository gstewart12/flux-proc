
library("units")
library("jsonlite")
library("glue")
library("tidyverse")
library("ecmwfr")

# Set token to local keychain
# Requires an ECMWF account 
# 1. Register at https://www.ecmwf.int/user/login
# 2. Get API token from CDS https://cds.climate.copernicus.eu/profile
# 3. Add token to keychain:
# wf_set_key(key = "")
#    OR input key interactively:
# wf_set_key()
# more info on setting up API here: https://cds.climate.copernicus.eu/how-to-api

area <- c(39.10, -75.80, 39.05, -75.75) # (N, W, S, E) includes all sites

era_vars <- list(
  instant = c(
    "10m_u_component_of_wind", 
    "10m_v_component_of_wind", 
    "2m_dewpoint_temperature", 
    "2m_temperature", 
    "surface_pressure",
    "friction_velocity",
    glue("volumetric_soil_water_layer_{c(1,2,3,4)}"),
    glue("soil_temperature_level_{c(1,2,3,4)}"),
    "boundary_layer_height"
  ),
  average = c(
    "mean_total_precipitation_rate", 
    "mean_surface_downward_long_wave_radiation_flux",
    "mean_surface_downward_short_wave_radiation_flux",
    "mean_surface_net_long_wave_radiation_flux",
    "mean_surface_net_short_wave_radiation_flux",
    "mean_surface_sensible_heat_flux",
    "mean_surface_latent_heat_flux"
  )
)


# Download one year of ERA5 data from CDS --------------------------------------

# Set ERA5 data directory
era_dir <- "data/ERA5"

year <- 2025

# Create data query/queries
req <- list(
  dataset_short_name = "reanalysis-era5-single-levels",
  product_type = "reanalysis",
  # variable = era_vars,
  date = paste0(year, "-01-01/", year, "-12-31"),
  time = paste0(str_pad(0:23, 2, pad = "0"), ":00"),
  data_format = "netcdf",
  download_format = "unarchived",
  area = area,
  target = paste0("era5-sl-", year, ".nc")
)

reqs <- map_depth(
  era_vars, 2, \(x) assign_in(req, "variable", x), .ragged = TRUE
)
reqs$instant <- map2(
  reqs$instant, 1:length(reqs$instant),
  \(x, i) modify_in(
    x, "target", 
    \(y) str_replace(y, "\\.", str_c("-", str_pad(i, 2, pad = "0"), "."))
  )
)
reqs$average <- map2(
  reqs$average, 1:length(reqs$average) + length(reqs$instant),
  \(x, i) modify_in(
    x, "target", 
    \(y) str_replace(y, "\\.", str_c("-", str_pad(i, 2, pad = "0"), "."))
  )
)

# Check size of request, split into multiple if needed
# req_size_limit <- 120000
# Handle instants and averages separately - cannot be mixed in a single request
# reqs <- list()
# for (i in 1:length(era_vars)) {
#   vars <- era_vars[[i]]
#   req_size <- prod(
#     length(req$time),
#     as.numeric(as.duration(interval(req$date)), "days"),
#     length(vars),
#     6 # correction factor for when request size suddenly grew 6x
#   )
#   req_max_vars <- req_size_limit %/% (req_size/length(vars))
#   n_reqs <- ceiling(length(vars)/req_max_vars)
#   req_vars <- split(vars, cut(seq_along(vars), n_reqs, labels = FALSE))
#   reqs[[names(era_vars[i])]] <- req_vars |>
#     unname() |>
#     map(\(x) assign_in(req, "variable", x))
# }
# reqs <- reqs |>
#   list_c() |>
#   imap(
#     \(x, i) modify_in(
#       x, "target", 
#       \(y) str_replace(y, "\\.", str_c("-", str_pad(i, 2, pad = "0"), "."))
#     )
#   )

# Split requests into two batches - ECMWF allows max 20 at once
batches <- list()
req_files <- list()

# First batch: instants
# Request the data - don't download yet, takes a while to process
batches$instant <- map(reqs$instant, \(x) wf_request(x, transfer = FALSE))

# Add request IDs to queries and save
reqs$instant <- map2(
  reqs$instant, batches$instant, \(x, y) append(x, list(url = y$get_url()))
)
req_files$instant <- reqs$instant |>
  map("target") |>
  map(\(x) str_replace(x, "\\.nc", "-req.json"))
walk2(
  reqs$instant, req_files$instant, \(x, y) write_json(x, file.path(era_dir, y))
)

# STOP, wait here until first batch is ready
# Check status: https://cds.climate.copernicus.eu/requests?tab=all
dminutes(9.5 * length(batches$instant)) # est. processing time (mins)

# Second batch: averages
# Request the data - don't download yet, takes a while to process
batches$average <- map(reqs$average, \(x) wf_request(x, transfer = FALSE))

# Add request IDs to queries and save
reqs$average <- map2(
  reqs$average, batches$average, \(x, y) append(x, list(url = y$get_url()))
)
req_files$average <- reqs$average |>
  map("target") |>
  map(\(x) str_replace(x, "\\.nc", "-req.json"))
walk2(
  reqs$average, req_files$average, \(x, y) write_json(x, file.path(era_dir, y))
)

# STOP, wait here until second batch is ready
# Check status: https://cds.climate.copernicus.eu/requests?tab=all
dminutes(9.5 * length(batches$average)) # est. processing time (mins)

# Download requested data ------------------------------------------------------

# Check status: https://cds.climate.copernicus.eu/requests?tab=all

reqs <- era_dir |>
  list.files(glue("{year}(-\\d{{1,2}})*-req"), full.names = TRUE) |>
  map(\(x) read_json(x, simplifyVector = TRUE))
walk(
  reqs,
  \(x) wf_transfer(
    url = x$url, 
    path = era_dir, 
    filename = x$target
  )
)

# Process ERA5 for use in tower gap filling ------------------------------------

# Read in data
nc_files <- list.files(era_dir, pattern = ".nc", full.names = TRUE)
nc_data <- nc_files |> 
  as.list() |>
  # Group files by year
  split(str_extract(nc_files, "(\\d{4})")) |>
  map_depth(2, quietly(stars::read_ncdf)) |>
  modify_depth(2, "result") |>
  map_depth(2, as_tibble) |>
  # Join multiple files from the same year
  map(\(x) reduce(x, full_join)) |> 
  # Recent data has two 'experiment versions' - need to remove one
  map(\(x) filter(x, if_any(any_of("expver"), \(x) x != 3))) |>
  map(\(x) select(x, -any_of("expver"))) |>
  # Join all years
  bind_rows() |>
  mutate(
    time = coalesce(time, valid_time),
    mtpr = coalesce(mtpr, avg_tprate),
    msdwlwrf = coalesce(msdwlwrf, avg_sdlwrf),
    msdwswrf = coalesce(msdwswrf, avg_sdswrf),
    msnlwrf = coalesce(msnlwrf, avg_snlwrf),
    msnswrf = coalesce(msnswrf, avg_snswrf),
    msshf = coalesce(msshf, avg_ishf),
    mslhf = coalesce(mslhf, avg_slhtf),
    .keep = "unused"
  ) |>
  drop_na() |>
  select(-longitude, -latitude) |>
  distinct(time, .keep_all = TRUE) |>
  arrange(time)

# Look at variable attributes
# nc_files |>
#   map(ncmeta::nc_atts) |>
#   map(\(x) select(x, -id)) |>
#   map(\(x) mutate(x, value = map_chr(value, as.character))) |>
#   map(pivot_wider) |>
#   bind_rows() |>
#   type_convert(col_types = cols()) |>
#   distinct(variable, .keep_all = TRUE) |>
#   print(n = 30)

# Upsample to flux tower time step
time_hh <- seq(
  min(nc_data$time) + minutes(15), 
  max(nc_data$time) + minutes(45), 
  by = "30 mins"
)

# ERA5 rate/flux parameters are averaged over the hour ending at 'time'
# - so the midpoint time for these variables is really time - minutes(30)
rate_vars <- c(
  "mtpr", "msdwlwrf", "msdwswrf", "msnlwrf", "msnswrf", "msshf", "mslhf"
)

nc_data_hh <- as.list(nc_data[0, ])
nc_data_hh[[1]] <- time_hh
for (i in 2:length(nc_data_hh)) {
  # Using Akima splines - smooth interpolation through given points
  if (names(nc_data_hh)[i] %in% rate_vars) {
    var_hh <- akima::aspline(nc_data$time - minutes(30), nc_data[[i]], time_hh)
  } else {
    var_hh <- akima::aspline(nc_data$time, nc_data[[i]], time_hh)
  }
  # Re-attach units
  nc_data_hh[[i]] <- as_units(var_hh$y, units(nc_data[[i]]))
}

era_data_hh <- nc_data_hh |>
  as_tibble() |>
  # Rename variables to match EddyPro/Fluxnet conventions
  rename(
    TIMESTAMP = time, TDEW = d2m, TA = t2m, LW_IN = msdwlwrf, SW_IN = msdwswrf, 
    PA = sp, H = msshf, LE = mslhf, BLH = blh, USTAR = zust
  ) |>
  rename_with(\(x) str_replace(x, "stl", "TS_"), starts_with("stl")) |>
  rename_with(\(x) str_replace(x, "swvl", "SWC_"), starts_with("swvl")) |>
  # Convert units
  mutate(
    # Timestamp in local (eastern) time; set to middle of avg. period
    TIMESTAMP = force_tz(with_tz(TIMESTAMP, "Etc/GMT+5"), "UTC"),
    # Temperatures from K to C
    across(c(TDEW, TA, starts_with("TS_")), \(x) set_units(x, degree_Celsius)),
    # Soil water content from m3 m-3 to percent
    across(c(starts_with("SWC_")), \(x) set_units(x, `%`)),
    # Precipitation from kg m-2 s-1 to mm
    P_RAIN = mtpr * set_units(1800, s) * set_units(1, "mm/kg/m^-2"),
    # Pressure from Pa to kPa
    PA = set_units(PA, kPa),
    # H and LE: reverse sign so positive fluxes are away from surface
    across(c(H, LE), `-`),
    .keep = "unused"
  ) |>
  # Calculate derived variables
  mutate(
    # Outgoing radiation components
    LW_OUT = LW_IN - msnlwrf,
    SW_OUT = SW_IN - msnswrf,
    # Set some negative values to 0 (consequence of upsampling with splines)
    across(c(SW_IN, SW_OUT), \(x) pmax(x, set_units(0, W/m^2))),
    P_RAIN = pmax(P_RAIN, set_units(0, mm)),
    # PPFD from SW_IN w/ conversion factor
    PPFD_IN = set_units(bigleaf::Rg.to.PPFD(drop_units(SW_IN)), µmol/m^2/s),
    # Net radiation
    NETRAD = (SW_IN - SW_OUT) + (LW_IN - LW_OUT),
    # Soil heat flux
    G = NETRAD - (H + LE),
    # Humidity variables
    e = bigleaf::Esat.slope(drop_units(TDEW))$Esat,
    VPD = set_units(bigleaf::e.to.VPD(e, drop_units(TA)) * 10, hPa),
    RH = set_units(bigleaf::e.to.rH(e, drop_units(TA)) * 100, `%`),
    SPECIFIC_HUMIDITY = set_units(bigleaf::e.to.q(e, drop_units(PA)), kg/kg),
    # Wind speed & direction
    WD = set_units(
      (270 - drop_units(atan2(v10, u10)) * 180/pi) %% 360, degrees
    ),
    # Scale wind speed to tower height (~3 m)
    WS = sqrt(u10^2 + v10^2) * (log(3) / log(10))
  ) |>
  select(-e, -u10, -v10, -msnlwrf, -msnswrf)

# Visual checks
era_data_hh |>
  ggplot(aes(TIMESTAMP, TA)) +
  geom_line()

# Write processed ERA5 data to files (yearly)
era_years <- unique(year(era_data_hh$TIMESTAMP))
year_seqs <- era_years |>
  map(\(x) c(glue("{x}-01-01 00:15"), glue("{x}-12-31 23:45"))) |>
  map(ymd_hm) |>
  map(\(x) seq(x[1], x[2], by = "30 mins")) |>
  map(\(x) as_tibble_col(x, "TIMESTAMP"))
era_years |>
  map(\(x) filter(era_data_hh, year(TIMESTAMP) == x)) |>
  # Make a full year sequence
  map2(year_seqs, \(x, y) left_join(y, x, by = join_by(TIMESTAMP))) |>
  walk2(era_years, \(x, y) write_csv(x, glue("output/ERA5/era5-sl-{y}-hh.csv")))


