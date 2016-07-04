#' Calculate Weight of Evidence
#'
#' Calculates the Weight of Evidence, namely \code{log((n_good/Population)/(n_bad/Population))}, while handling special cases (either no goods, no bads or both). If all the buckets are special, the value defaults to \code{abs(extremum.woe)} if there are only goods, \code{-abs(extremum.woe)} if there are only goods and 0 if there are neither. If only a few buckets have no goods or bads, then the values default to \code{na.factor} times the smallest (if no goods) or highest present WoE(if no bads).
#'
#' @param n_good Number of negative examples ("goods")
#' @param n_good Number of positive examples ("bads")
#' @param na.factor Multiplicative facor in case some cases are NAs
#' @param extremum.woe WoE in case all observations are NA
#' @export
woe <- function(n_good, n_bad, na.factor = 2, extremum.woe = 5){
  w <- ifelse(n_good == 0 | n_bad == 0 | is.na(n_good) | is.na(n_bad),
         NA,
         log(n_good/n_bad))

  if(any(is.na(w))){
    if(all(is.na(w))){
      message('All buckets have observations of a single class. Replacing WoE with +5 if there are only goods and -5 if there are only bads.')
      w <- ifelse(n_bad == 0,
                  ifelse(n_good == 0,
                         0,
                         abs(extremum.woe)),
                  -abs(extremum.woe))
    } else{
      message('Some buckets have observations of a single class. Replacing WoE with twice the maximum / minimum WoE among other buckets.')
      w <- ifelse(is.na(w),
                  ifelse(n_bad == 0,
                         na.factor*max(w, na.rm = T),
                         na.factor*min(w, na.rm = T)),
                  w)
    }
  }
  return(w)
}

# Local version of %>% in case it is not already imported and so that I don't need to use dplyr::%>%
#' @export
`%>%` <- dplyr::`%>%`

# Convert integer columns in a data.frame to double
intcols2double <- function(df){
  cl <- class(df)
  ix <- sapply(df, is.integer)
  df[,ix] <- sapply(df[,ix], as.double)
  df
}

# cut has a bug with frontier values sometimes
cut_ <- function(x, breaks){
  breaks <- sort(unique(breaks))
  mins <- breaks[-length(breaks)]
  maxs <- breaks[-1]
  out <- sapply(x, function(y){
    cond <- which(y > mins & y <= maxs)
    ifelse(length(cond) == 0,
           NA,
           cond)
  })
  factor(out,
         levels = 1:(length(breaks) - 1),
         labels = sprintf('(%.2f,%.2f]', mins, maxs))
}

qcut <- function(x, g=10){
  if(!is.numeric(g) || g <= 0) stop("g should be a positive integer.")
  if(is.factor(x)){
    warning("Converting factor to numeric.")
    x <- as.numeric(x)
  }
  percentages <- seq(1/g, 1-(1/g), l=g-1)
  qq <- quantile(x, percentages)
  y <- cut_(x,
            breaks=c(-Inf,
                     unique(qq),
                     Inf))
  list(
    x = y,
    percentages = c(round(percentages, 3), 1)
  )
}

# Based on included cut()
qcut2 <- function(x, g=10){
  if(!is.numeric(g) || g <= 0) stop("g should be a positive integer.")
  if(is.factor(x)){
    warning("Converting factor to numeric.")
    x <- as.numeric(x)
  }
  percentages <- seq(1/g, 1-(1/g), l=g-1)
  qq <- quantile(x, percentages)
  y <- cut(x,
           breaks=c(-Inf,
                    unique(qq),
                    Inf),
           include.lowest = TRUE,
           right = TRUE)
  list(
    x = y,
    percentages = c(round(percentages, 3), 1)
  )
}

# Using weighted quantiles
qcut3 <- function(x, wt, g=10){
  if(!is.numeric(g) || g <= 0) stop("g should be a positive integer.")
  if(is.factor(x)){
    warning("Converting factor to numeric.")
    x <- as.numeric(x)
  }
  percentages <- seq(1/g, 1-(1/g), l=g-1)
  qq <- Hmisc::wtd.quantile(x, weights = wt, probs = percentages)
  y <- cut(x,
           breaks=c(-Inf,
                    unique(qq),
                    Inf),
           include.lowest = TRUE,
           right = TRUE)
  list(
    x = y,
    percentages = c(round(percentages, 3), 1)
  )
}
