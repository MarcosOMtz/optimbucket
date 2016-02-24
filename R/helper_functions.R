
require(dplyr)
require(tidyr)
require(ggplot2)

qcut <- function(x, g=10){
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
           right=TRUE)
  list(
    x = y,
    percentages = c(round(percentages, 3), 1)
  )
}
