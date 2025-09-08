
# Find all relevant EddyPro output
get_eddypro_files_year <- function(dir, year, type = "fluxnet") {
  
  # Finds EddyPro output files covering a given year
  
  year <- as.character(year)
  
  output_dirs <- stringr::str_subset(list.dirs(dir), "output$")
  
  # output_dates <- output_dirs |>
  #   purrr::map(\(x) list.files(x, pattern = "_metadata_", full.names = TRUE)) |>
  #   purrr::map(\(x) paste("cut -f2 -d,", x)) |>
  #   purrr::map(\(x) scan(pipe(x), what = character(), skip = 1, quiet = TRUE))
  
  # Check range of dates in each EddyPro metadata file
  output_dates <- output_dirs |>
    purrr::map(\(x) list.files(x, pattern = "_metadata_", full.names = TRUE)) |>
    purrr::map(
      \(x) readr::read_csv(x, col_types = readr::cols_only(date = "c"))
    ) |>
    purrr::map(\(x) if (nrow(x)) dplyr::pull(x, date) else c(NA, NA))
  
  starts <- purrr::map_chr(output_dates, dplyr::first)
  ends <- purrr::map_chr(output_dates, dplyr::last)
  
  # Find files overlapping selected year
  is_year <- stringr::str_detect(starts, year) | stringr::str_detect(ends, year)
  year_dirs <- na.omit(output_dirs[is_year])
  
  files <- purrr::map_chr(year_dirs, \(x) list.files(x, pattern = type))
  file.path(year_dirs, files)
}

# get_debias_window <- function(x, y, target_ind, na_frac) {
#   stopifnot(length(x) == length(y))
#   w_start <- target_ind[1]
#   w_end <- target_ind[length(target_ind)]
#   n_hh <- w_end - w_start + 1
#   while (TRUE) {
#     x_w <- x[w_start:w_end]
#     y_w <- y[w_start:w_end]
#     # Stop if window has enough valid data
#     if (sum(complete.cases(x_w, y_w)) / n_hh >= na_frac) {
#       break
#     }
#     # Increase window by ±1 day
#     if (w_start > 1) w_start <- w_start - 48
#     if (w_end < length(x)) w_end <- w_end + 48
#     n_hh <- w_end - w_start + 1
#     if (n_hh == length(x)) break
#   }
#   w_start:w_end
# }

get_debias_window <- function(x, y, target_ind, na_frac) {
  stopifnot(length(x) == length(y))
  window_start <- target_ind[1]
  window_end <- target_ind[length(target_ind)]
  max_expand <- max(window_start - 1, length(y) - window_end)

  # Precompute valid indices
  is_valid <- complete.cases(x, y)
  valid_cumsum <- cumsum(is_valid)
  
  # Function to get valid count in a window [i, j]
  count_valid <- function(i, j) {
    if (i == 1) return(valid_cumsum[j])
    return(valid_cumsum[j] - valid_cumsum[i - 1])
  }
  
  # Binary search for shortest expansion satisfying max percent NA
  best_window <- NULL
  for (e in seq(0, max_expand, by = 48)) {
    i <- max(1, window_start - e)
    j <- min(length(y), window_end + e)
    valid_count <- count_valid(i, j)
    if (valid_count / (j - i + 1) >= na_frac) {
      best_window <- i:j
      break
    }
  }
  best_window
}

# Function to compute lagged correlation and detect phase shift
get_best_lag <- function(x, y, lag.max = 48, lag.default = 0) {
  
  if (all(is.na(x + y))) {
    return(lag.default)
  }
  
  ccf <- ccf(x, y, na.action = na.pass, lag.max = lag.max, plot = FALSE)
  
  # return k, where x[t+k] ~ y[t]
  ccf$lag[, , 1][which.max(ccf$acf[, , 1])]
}

apply_lag <- function(x, k) {
  if (k > 0) {
    # if k > 0, x lags y (changes in x follow changes in y)
    dplyr::lag(x, n = k)
  } else {
    # if k < 0, x leads y (changes in y follow changes in x)
    dplyr::lead(x, n = -k)
  }
}

# [function purpose]
cov_wd <- function(x, y, deg = 0:359) {
  map_dbl(deg, \(d) cov((x + d) %% 360, y, use = "pairwise.complete.obs"))
}

get_fit_stats <- function(x, y, fit) {
  n_hh <- length(x)
  
  if (inherits(fit, "lm")) {
    list(
      frac_hh = length(fitted(fit))/n_hh,
      slope = if (is.na(coef(fit)[2])) coef(fit)[1] else coef(fit)[2],
      int = if (is.na(coef(fit)[2])) 0 else coef(fit)[1],
      r2 = summary(fit)$r.squared,
      rmse = sqrt(mean(fit$residuals^2)),
      # var_ratio = var(fit$model[, 1])/var(fit$model[, 2]),
      var_ratio = var(fit$model[, 1])/var(fitted(fit)),
      n_used = 0
    )
  } else if (inherits(fit, "numeric")) {
    int <- which.max(fit) - 1
    list(
      frac_hh = sum(complete.cases(x, y))/n_hh,
      slope = 1,
      int = int,
      r2 = cor((x + int) %% 360, y, use = "pairwise.complete.obs")^2,
      rmse = NA,
      var_ratio = NA,
      n_used = 0
    )
  }
}

