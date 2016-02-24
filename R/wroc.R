
# PRUEBAS TRAMEADOS

require(dplyr)
require(tidyr)
require(ggplot2)

wroc <- function(x, ...) UseMethod("wroc")

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

reset.buckets <- function(x){
  x$info$bucket <- 0:(nrow(x$info)-1)
  x
}

performance.wroc.info_ <- function(ds, trend = c('upper','lower')){

  if(!(trend[1] %in% c('upper','lower'))){
    warning('Not a valid trend. Only "upper" and "lower" are allowed. Defaulting to "upper".')
    trend <- 'upper'
  }
  sign_trend <- ifelse(trend[1]=='upper', 1, -1)

  tr <- ds %>%
    mutate(dx = d_ac_good - dplyr::lag(d_ac_good),
           dx_skip = dplyr::lead(d_ac_good) - dplyr::lag(d_ac_good),
           bases_gini_raw = pmax(sign_trend*(d_ac_bad - d_ac_good), 0),
           bases_gini = bases_gini_raw + lag(bases_gini_raw),
           bases_gini_skip = lead(bases_gini_raw) + lag(bases_gini_raw),
           area_gini = bases_gini*dx/2,
           area_gini_skip = bases_gini_skip*dx_skip/2,
           delta_area = area_gini + lead(area_gini)- area_gini_skip)

  tr
}

choose.trend_ <- function(x){
  n_below <- sum(x$info$d_ac_bad < x$info$d_ac_good)
  trend <- ifelse(n_below >= nrow(x$info)/2, 'lower', 'upper')
  trend
}

wroc.default <- function(predictions, labels, ngroups=50, level.bad=1, col.bad=1,
                         special.values = NULL){
  if(!is.numeric(predictions)){
    warning("Predictions should be numeric. Coercing to numeric.")
    predictions <- as.numeric(predictions)
  }

  out <- list()
  out$call <- match.call()
  class(out) <- 'wroc'

  if(!is.null(special.values)){
    special.values <- unique(special.values)
    special_ix <- sapply(predictions, function(i){
      (i %in% special.values)
    })
    special_buckets <- sapply(predictions[special_ix], function(i){
      -which(i == special.values)
    })
    special_predictions <- predictions[special_ix]
    special_labels <- labels[special_ix]
  } else {
    special_ix <- rep(F, length(predictions))
    special_buckets <- NULL
  }


  buckets <- rep(NA, length(special_ix))
  buckets[special_ix] <- special_buckets
  if(is.na(ngroups)){
    warning('Using exact ROC curve. This may be very slow for continuous variables! Try using a smaller number for ngroups.')
    buckets[!special_ix] <- 1:sum(!special_ix)
  } else{
    # buckets <- as.numeric(cut2(predictions, g = ngroups))
    buckets[!special_ix] <- as.numeric(qcut(predictions[!special_ix], g = ngroups)$x)
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


  ix <- (out$info$bucket %in% special_buckets)
  out$special <- out$info[ix,]
  out$info <- out$info[!ix,]
  out$info$upper_limit[1] <- -Inf
  out$info$upper_limit[nrow(out$info)] <- Inf
  out$info$lower_limit <- lag(out$info$upper_limit)
  out$info$lower_limit[1] <- -Inf


  out$ngroups <- nrow(out$info) - 1
  out$nspecial <- nrow(out$special)

  out
}

plot.wroc <- function(x,
                      type = c('accum','roc','trend','woe'),
                      include.special = TRUE){
  require(ggplot2)

  if(type[1] == 'accum'){
    p <- x$info %>%
      ggplot(aes(d_ac_good, d_ac_bad, color=bucket)) +
      geom_rect(xmin=0,xmax=1,ymin=0,ymax=1,fill=NA, color='black') +
      geom_segment(aes(x=0,xend=1,y=0,yend=1), color='black', linetype='dashed') +
      geom_line(size=1) +
      geom_point() +
      geom_point(size=1.5, shape=1, color='black') +
      scale_color_gradientn(colors = c('blue','green','yellow','red'))
  } else if(type[1] == 'roc'){
    p <- x$info %>%
      ggplot(aes(1-d_ac_good, 1-d_ac_bad, color=bucket)) +
      geom_rect(xmin=0,xmax=1,ymin=0,ymax=1,fill=NA, color='black') +
      geom_segment(aes(x=0,xend=1,y=0,yend=1), color='black', linetype='dashed') +
      geom_line(size=1) +
      geom_point() +
      geom_point(size=1.5, shape=1, color='black') +
      scale_color_gradientn(colors = c('blue','green','yellow','red')) +
      labs(x='fpr',y='tpr')
  } else if(type[1] == 'trend'){
    ds <- x$info[-1,]
    if(include.special) ds <- rbind(x$special, ds)
    ds <- ds %>%
      mutate(barcol = ifelse(bucket < 0, 'red', 'darkgrey'))
    brks <- 1:nrow(ds)
    labls <- sprintf('B%d: (%.2f, %.2f]',
                     ds$bucket,
                     ds$lower_limit,
                     ds$upper_limit)
    if(include.special) labls[1:nrow(x$special)] <- gsub('\\(','[',labls[1:nrow(x$special)])
    labls[length(labls)] <- gsub(']',')',labls[length(labls)])
    p <- ds %>%
      mutate(i = row_number(),
             norm_population = population*max(p_bad)/max(population)) %>%
      ggplot(aes(i, p_bad)) +
      geom_bar(aes(y=norm_population, fill=barcol), stat='identity') +
      geom_point() +
      geom_line() +
      geom_text(aes(y = 0, label=sprintf('%.2f %%',d_population)),
                size = 2, vjust=1)+#angle=90, hjust = -0.5) +
      scale_fill_identity() +
      scale_x_continuous(breaks=brks,labels=labls) +
      theme(axis.text.x = element_text(angle=90)) +
      labs(x = 'Bucket',
           y = 'Default rate')

  } else if(type[1] == 'woe'){
    ds <- x$info[-1,]
    if(include.special) ds <- rbind(x$special, ds)
    ds <- ds %>%
      mutate(barcol = ifelse(bucket < 0, 'red', 'darkgrey'))
    brks <- 1:nrow(ds)
    labls <- sprintf('B%d: (%.2f, %.2f]',
                     ds$bucket,
                     ds$lower_limit,
                     ds$upper_limit)
    if(include.special) labls[1:nrow(x$special)] <- gsub('\\(','[',labls[1:nrow(x$special)])
    labls[length(labls)] <- gsub(']',')',labls[length(labls)])
    p <- ds %>%
      mutate(i = row_number(),
             norm_population = population*max(woe)/max(population)) %>%
      ggplot(aes(i, woe)) +
      geom_bar(aes(y=norm_population, fill=barcol), stat='identity') +
      geom_point() +
      geom_line() +
      geom_text(aes(y = 0, label=sprintf('%.2f %%',d_population)),
                size = 2, vjust=1)+#angle=90, hjust = -0.5) +
      scale_fill_identity() +
      scale_x_continuous(breaks=brks,labels=labls) +
      theme(axis.text.x = element_text(angle=90)) +
      labs(x = 'Bucket',
           y = 'Weight of Evidence')
  }
 p
}
optimize.wroc <- function(x, trend = c('auto','upper','lower')){
  if('optimal.wroc' %in% class(x)){
    warning('This is an already optimal ROC curve. Returning the input ROC curve. ')
    return(x)
  } else if('subset.wroc' %in% class(x)){
    warning('This ROC curve has been manually tampered with. Resetting bucket names.')
    x <- reset.buckets(x)
  }
  ds <- x$info
  if(!(trend[1] %in% c('auto','upper','lower'))){
    warning('Invalid trend. Only "upper", "lower" or "auto" are valid. Defaulting to "auto".')
    trend <- 'auto'
  }
  if(trend[1] == 'auto'){
    trend <- choose.trend_(x)
  }

  which.fun <- ifelse(trend[1] == 'upper', which.max, which.min)

  jumps <- sapply(1:(nrow(ds)-1), function(i){
    len <- nrow(ds)
    dx <- ds$d_ac_good[(i+1):len] - ds$d_ac_good[i]
    dy <- ds$d_ac_bad[(i+1):len] - ds$d_ac_bad[i]
    j <- ifelse(all(dy == 0) || all(dx == 0), length(dy),
                which.fun(ifelse(dx == 0, Inf, dy/dx)))
    j
  })

  ix <- rep(1, length(jumps)+1)
  i <- 1
  while(i <= length(jumps)){
    if(jumps[i] > 1){
      ix[i+1] <- jumps[i]
      ix[(i+2):(i+jumps[i])] <- 0
      i <- i + jumps[i]
    } else{
      i <- i + 1
    }
  }

  aux <- cumsum(ix) - 1
  ixx <- unique(aux)
  buckets_to_remove <- which(!(1:x$ngroups %in% ixx))
  out <- subset(x, buckets = buckets_to_remove)

  out$removed.buckets <- c(x$removed.buckets, buckets_to_remove)
  out$trend <- trend
  out$call.optimize <- match.call()

  class(out) <- c('optimal.wroc', 'wroc')
  out
}

performance.wroc <- function(x){
  ds <- x$info

  tr <- ds %>%
    mutate(dx = (d_ac_good - dplyr::lag(d_ac_good)),
           bases_auc = (d_ac_bad + dplyr::lag(d_ac_bad)),
           bases_gini_up_raw = pmax(d_ac_bad - d_ac_good, 0),
           bases_gini_up = bases_gini_up_raw + lag(bases_gini_up_raw),
           bases_gini_down_raw = pmax(d_ac_good - d_ac_bad, 0),
           bases_gini_down = bases_gini_down_raw + lag(bases_gini_down_raw),
           area_auc = bases_auc*dx/2,
           area_gini_up = bases_gini_up*dx/2,
           area_gini_down = bases_gini_down*dx/2,
           ks_prospect = abs(d_good - d_bad),
           iv_contribution = (d_good - d_bad)*woe) %>%
    .[-1,]

  out <- list()
  out$auc <- sum(tr$area_auc)
  out$gini_two_way <- abs(2*out$auc - 1)
  out$gini_up <- 2*sum(tr$area_gini_up)
  out$gini_down <- 2*sum(tr$area_gini_down)
  out$gini <- max(out$gini_up, out$gini_down)
  out$ks <- max(tr$ks_prospect)
  out$iv <- sum(tr$iv_contribution)

  class(out) <- c('wroc.performance')
  out
}

subset.wroc <- function(x, buckets = NULL, ...){
  if(is.null(buckets)) return(x)
  if(any(buckets < 0)){
    stop('Cannot paste special valued buckets. Join the levels by hand before calling wroc.')
  } else {
    ds <- x$info[-1,]
  }

  if(0 %in% buckets){
    buckets <- buckets[buckets > 0]
    warning('Buckets must be strictly positive integers.')
  }
  if(max(ds$bucket) %in% buckets){
    stop("Cannot remove the last bucket.")
  }

  buckets <- sort(unique(buckets))
  ix <- sort(sapply(buckets, function(y) which(y == ds$bucket)))
  # This must NOT be vectorized or else it doesn't work
  for(ii in rev(ix)){
    ds$bucket[ii] <- ds$bucket[ii+1]
  }

  out <- wroc(predictions=ds$bucket,
              labels=cbind(ds$n_bad, ds$n_good),
              ngroups = nrow(ds))

  out$special <- x$special
  out$nspecial <- x$nspecial
  out$info$bucket <- c(0, ds$bucket[-ix])
  out$info$lower_limit <- c(-Inf, ds$lower_limit[-(ix+1)])
  out$info$upper_limit <- c(-Inf, ds$upper_limit[-(ix)])
  out$ngroups <- nrow(out$info) - 1
  out$call.subset <- match.call()
  out$pasted_buckets <- buckets
  class(out) <- c('subset.wroc', class(out))
  out
}


analyze.wroc <- function(x,
                         nbuckets = 5,
                         trend = c('auto','upper','lower')){

  if(is.null(nbuckets) || nbuckets < 0 || !is.numeric(nbuckets)){
    message('analyze.wroc: Invalid number of buckets. Optimizing instead.')
    return(optimize.wroc(x, trend))
  }

  if(!(trend[1] %in% c('auto','upper','lower'))){
    warning('Invalid trend. Only "upper", "lower" or "auto" are valid. Defaulting to "auto".')
    trend <- 'auto'
  }
  if(trend[1] == 'auto'){
    trend <- choose.trend_(x)
  }

  if(x$ngroups < nbuckets){
    stop("Cannot remove more buckets than there are!")
  }

  ds <- x$info

  tr <- performance.wroc.info_(x$info, trend = trend)


  buckets_to_paste <- numeric(x$ngroups - nbuckets)
  gini_loss <- numeric(x$ngroups - nbuckets)
  gini <- numeric(x$ngroups - nbuckets)
  for(k in 1:(x$ngroups-nbuckets)){
    j <- which.min(tr$delta_area)
    buckets_to_paste[k] <- tr$bucket[j]
    tr <- tr[-j,]
    if(j == 2){
      aux <- tr[1:3,]
      tr[2,] <- performance.wroc.info_(aux, trend = trend)[2,]
    } else if(j > 2 && j < nrow(ds)){
      aux <- tr[(j-2):(j+1),]
      tr[(j-1):j,] <- performance.wroc.info_(aux, trend = trend)[2:3,]
    }
    gini[k] <- 2*sum(tr$area_gini, na.rm = T)
    if(k == 1){
      if(trend[1] == 'upper'){
        gini_at_start <- performance.wroc(x)$gini_up
      } else {
        gini_at_start <- performance.wroc(x)$gini_down
      }
      gini_loss[k] <- gini_at_start - gini[k]
    } else {
      gini_loss[k] <- gini[k-1] - gini[k]
    }
  }

  out <- subset(x, buckets = buckets_to_paste)
  out$trend <- trend
  out$pasted_buckets <- buckets_to_paste
  out$gini_at_start <- gini_at_start
  out$gini <- gini
  out$gini_loss <- gini_loss
  out$class <- c('analyze.wroc', class(out))
  out
}

wroc.formula <- function(formula, data, ngroups = 50, level.bad=1,
                         special.values=NULL){
  if(is.null(special.values)){
    special.values <- list()
  } else if((!(is.list(special.values)) && !is.numeric(special.values))
            ||
            (is.list(special.values) && is.null(names(special.values)))){
    stop('special.values must be either a named list with a vector of special values for each variable or a single vector of special values common to all the variables.')
  }

  reuse_special <- (is.numeric(special.values) || is.null(special.values))

  ds <- model.frame(formula, data)
  y <- ds[[1]]
  ds <- ds[-1]

  out <- lapply(1:ncol(ds), function(j){
    if(reuse_special){
      spvals <- special.values
    } else {
      spvals <- special.values[[names(ds)[j]]]
    }
    wroc(ds[[j]], y, ngroups = ngroups, level.bad = level.bad,
         special.values = spvals)
  })
  names(out) <- names(ds)
  class(out) <- 'wroc.conglomerate'
  out
}

wrs <- wroc(y ~ x + z, d, ngroups = 20, level.bad = 1)

#compactify.wroc






