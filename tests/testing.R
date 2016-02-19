wroc.default <- function(predictions, labels, ngroups=50, level.bad=1, col.bad=1){
  if(!is.numeric(predictions)){
    warning("Predictions should be numeric. Coercing to numeric.")
    predictions <- as.numeric(predictions)
  }

  out <- list()
  out$call <- match.call()
  class(out) <- 'wroc'

  if(is.na(ngroups)){
    buckets <- 1:length(predictions)
  } else{
    # buckets <- as.numeric(cut2(predictions, g = ngroups))
    buckets <- as.numeric(qcut(predictions, g = ngroups)$x)
  }

  if((is.vector(labels) && is.numeric(labels)) || is.factor(labels)){
    out$info <- data.frame(bucket = buckets,
                           x = predictions,
                           y = labels) %>%
      group_by(bucket) %>%
      summarise(population = n(),
                lower_limit = max(x),
                upper_limit = max(x),
                n_bad = sum(y == level.bad))
  } else if((is.matrix(labels) | is.data.frame(labels)) && ncol(labels) == 2){
    ix <- (apply(labels, 1, sum) != 0)
    out$info <- data.frame(bucket = buckets,
                           x = predictions,
                           n_bad = labels[,col.bad],
                           n_good = labels[,which(1:2 != col.bad)]) %>%
      filter(ix) %>%
      group_by(bucket) %>%
      summarise(population = sum(n_bad) + sum(n_good),
                lower_limit = max(x),
                upper_limit = max(x),
                n_bad = sum(n_bad))
  } else{
    stop("The labels must be either a vector of classes (numeric or factor) or a two-column matrix with the counts for each class.")
  }

  out$info <- out$info %>%
    mutate(n_good = population - n_bad,
           p_bad = n_bad / population,
           p_good = 1 - p_bad,
           ac_population = cumsum(population),
           ac_bad = cumsum(n_bad),
           ac_good = cumsum(n_good),
           tot_population = sum(population),
           tot_bad = sum(n_bad),
           tot_good = sum(n_good),
           d_population = population/tot_population,
           d_bad = n_bad/tot_bad,
           d_good = n_good/tot_good,
           d_ac_population = ac_population/tot_population,
           d_ac_bad = ac_bad/tot_bad,
           d_ac_good = ac_good/tot_good,
           woe = ifelse(p_bad == 0 | p_good == 0, NA, log(p_good/p_bad))) %>%
    rbind(0, .)

  if(any(is.na(out$info$woe))){
    warning('Some buckets have observations of a single class. Replacing WoE with twice the maximum / minimum WoE among other buckets.')
    out$info <- out$info %>%
      mutate(woe = ifelse(is.na(woe),
                          ifelse(p_bad == 0,
                                 2*max(woe, na.rm = T),
                                 2*min(woe, na.rm = T)),
                          woe))
  }


  out$info$upper_limit[1] <- -Inf
  out$info$upper_limit[nrow(out$info)] <- Inf
  out$info$lower_limit <- lag(out$info$upper_limit)
  out$info$lower_limit[1] <- -Inf

  out$ngroups <- nrow(out$info) - 1

  out
}
