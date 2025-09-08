# Flexible method for flagging time series spikes
# After testing many methods, this one seems most efficient across all variables
# This is a re-implementation of the method in NEON's processing pipeline
# With a couple of changes e.g. using the Qn estimator instead of MAD
flag_spikes <- function(x, z = 5, lim_flag = NULL) {
  
  if (!is.null(lim_flag)) {
    x <- dplyr::if_else(lim_flag == 2, NA_real_, x)
  }
  
  # Calculate x(t) - 1/2*[x(t-1) + x(t+1)]
  y <- x - (lag(x) + lead(x))/2
  
  md <- median(y, na.rm = TRUE)
  # scale <- mad(y, na.rm = TRUE)
  scale <- robustbase::Qn(y, na.rm = TRUE)
  
  # Helper function for calculating excess kurtosis
  kurtosis <- function(x) {
    x <- na.omit(x)
    n <- length(x)
    # Pearson's definition of excess kurtosis
    g2 <- n * sum((x - mean(x))^4) / (sum((x - mean(x))^2)^2) - 3
    # Sample size adjustment
    (g2 + 3) * (1 - 1 / n)^2 - 3
  }
  
  # Calculate threshold
  kur <- kurtosis(y)
  a <- -0.2775258
  b <- 0.1720364
  z_adj <- exp(a - (b * log(kur)))/exp(a - (b * log(3)))
  
  z_adj <- 1 # get rid of this to actually use adjustment
  z <- z * z_adj
  
  # Check y against threshold
  diff <- y - md
  thr <- z * scale
  
  zdiff <- dplyr::case_when(
    diff > +thr ~ +1,
    diff < -thr ~ -1,
    TRUE ~ 0
  )

  # Classify outliers
  lo <- lag(zdiff) == -1 & zdiff == +1 & lead(zdiff) == -1
  hi <- lag(zdiff) == +1 & zdiff == -1 & lead(zdiff) == +1
  
  flag <- dplyr::if_else(hi | lo, 2L, 0L)
  tidyr::replace_na(flag, 0)
}

