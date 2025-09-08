
library("units")
library("jsonlite")
library("glue")
library("tidyverse")
library("ecmwfr")

area <- "39.10/-75.80/39.05/-75.75" # these bounds (N/W/S/E) include all sites
era_vars <- c(
  "10m_u_component_of_wind", 
  "10m_v_component_of_wind", 
  "2m_dewpoint_temperature", 
  "2m_temperature", 
  "boundary_layer_height",
  "mean_total_precipitation_rate", 
  "mean_surface_downward_long_wave_radiation_flux",
  "mean_surface_downward_short_wave_radiation_flux",
  "mean_surface_net_long_wave_radiation_flux",
  "mean_surface_net_short_wave_radiation_flux",
  "surface_pressure",
  "mean_surface_sensible_heat_flux",
  "mean_surface_latent_heat_flux",
  "friction_velocity",
  glue("volumetric_soil_water_layer_{c(1,2,3,4)}"),
  glue("soil_temperature_level_{c(1,2,3,4)}")
)


# Download one year of ERA5 data from CDS --------------------------------------

# Set ERA5 data directory
era_dir <- "data/ERA5"

# Set token to local keychain
wf_set_key(service = "cds") # open a browser window, enter key interactively
# Or enter credentials directly:
# wf_set_key(
#   user = "", 
#   key = "", 
#   service = "cds"
# )

year <- 2024

# Create data query/queries
req <- list(
  dataset_short_name = "reanalysis-era5-single-levels",
  product_type = "reanalysis",
  format = "netcdf",
  variable = era_vars,
  time = c(
    "00:00", "01:00", "02:00", "03:00", "04:00", "05:00", "06:00", "07:00", 
    "08:00", "09:00", "10:00", "11:00", "12:00", "13:00", "14:00", "15:00", 
    "16:00", "17:00", "18:00", "19:00", "20:00", "21:00", "22:00", "23:00"
  ),
  date = paste0(year, "-01-01/", year, "-12-31"),
  area = area,
  target = paste0("era5-sl-", year, ".nc")
)

# Check size of request, split into multiple if needed
req_size_limit <- 60000 # note this was lowered from 120000 to 60000
req_size <- prod(
  length(req$time), 
  as.numeric(as.duration(interval(req$date)), "days"),
  length(era_vars)
)
if (req_size > req_size_limit) {
  req_max_vars <- req_size_limit %/% (req_size/length(era_vars))
  n_reqs <- length(era_vars) %/% req_max_vars + 1
  req_vars <- split(era_vars, cut(seq_along(era_vars), n_reqs, labels = FALSE))
  reqs <- req_vars |> 
    map(\(x) assign_in(req, "variable", x)) |> 
    imap(
      \(x, i) modify_in(
        x, "target", \(y) str_replace(y, "\\.", str_c("-", i, "."))
      )
    )
} else {
  reqs <- list(req)
}

# Request the data - don't download yet, may take a while to process
files <- map(reqs, \(x) wf_request(x, user = "23795", transfer = FALSE))

# Add request IDs to queries and save
reqs <- map2(reqs, files, \(x, y) append(x, list(url = y$get_url())))
req_files <- reqs |>
  map("target") |>
  map(\(x) str_replace(x, "\\.nc", "-req.json"))
walk2(reqs, req_files, \(x, y) write_json(x, file.path(era_dir, y)))

dseconds(0.038) * req_size # Est. processing time (from 10 reps)

# STOP HERE

# Download requested data ------------------------------------------------------

# Check status: https://cds.climate.copernicus.eu/cdsapp#!/yourrequests

reqs <- map(
  list.files(era_dir, glue("{year}(-\\d)*-req"), full.names = TRUE), 
  \(x) read_json(x, simplifyVector = TRUE)
)
walk(
  reqs,
  \(x) wf_transfer(
    url = x$url, 
    user = "23795",
    path = era_dir, 
    filename = x$target, 
    service = "cds"
  )
)

# Process ERA5 for use in tower gap filling ------------------------------------

# Read in data
nc_files <- list.files(era_dir, pattern = ".nc", full.names = TRUE)
nc_data <- nc_files |> 
  as.list() |>
  # Group files by year
  split(str_sub(nc_files, -9, -6)) |>
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
  drop_na() |>
  select(-longitude, -latitude) |>
  distinct(time, .keep_all = TRUE) |>
  arrange(time)

# Look at variable attributes
nc_files |>
  map(ncmeta::nc_atts) |>
  map(\(x) select(x, -id)) |>
  map(\(x) mutate(x, value = map_chr(value, as.character))) |>
  map(pivot_wider) |>
  bind_rows() |>
  type_convert(col_types = cols()) |>
  distinct(variable, .keep_all = TRUE) |>
  print(n = 30)

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
  ggplot(aes(TIMESTAMP, P_RAIN)) +
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


