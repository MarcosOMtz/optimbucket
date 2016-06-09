
# Local helper functions

#' Reset Bucket Numbers of \code{wroc} Objects
#'
#' Resets the bucket numbers of an optimized or subsetted ROC curve in order to
#' have consecutive bucket numbers.
#'
#' @param x An object of type \code{wroc}.
#' @export
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
    dplyr::mutate(dx = d_ac_good - dplyr::lag(d_ac_good),
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

#' Create Flexible ROC Curves
#'
#' \code{wroc} builds ROC curves based on either a raw variable or a summarized
#' dataset. \code{wroc} objects can be manipulated in a variety of ways,
#' including manual and Optimal Gini bucketing.
#'
#' @param predictions A numeric vector or a matrix with one column with the
#'   predictions of a score or the value of a variable whose ROC curve is to be
#'   created.
#' @param labels Either a factor (or character vector) with the true labels of
#'   the response variable or a 2-column matrix containing counts for each class
#'   for each distinct value in \code{predictions}.
#' @param ngroups Precision to use for the ROC curve. A NULL value means "use as
#'   many values as there are unique values in \code{predictions}".
#' @param level.bad The level or value in \code{labels} corresponding to the
#'   positive class.
#' @param col.bad If \code{labels} is a matrix, the column containing the counts
#'   for the positive class.
#' @param special.values A vector containing values to be treated separately
#'   (NAs should be replaced by one such value).
#' @param formula Standard formula to specify the response (on the left side)
#'   and the variables whose ROC curves are to be created.
#' @param data A \code{data.frame} with the information.
#' @return An object of class \code{wroc} (default method) or a \code{wroc.list}
#'   of such objects (formula method and combine method), which can be then be
#'   manipulated in bulk.
#' @export
wroc <- function(x, ...) UseMethod("wroc")

#summarize_predictions_ <- function(predictions, labels, ngroups)

#' @describeIn wroc Method for raw, unsummarized data
#' @export
wroc.default <- function(predictions, labels, ngroups=NULL, level.bad=1, col.bad=1,
                         special.values = NULL){

  # Check types
  if(!is.numeric(predictions)){
    warning("Predictions should be numeric. Coercing to numeric.")
    predictions <- as.numeric(predictions)
  }
  if((is.vector(labels) && (is.numeric(labels) || is.character(labels))) ||
     is.factor(labels)){
    info_type <- 'long'
  } else if((is.matrix(labels) || is.data.frame(labels)) && ncol(labels) == 2){
    info_type <- 'summarized'
  } else{
    stop("The labels must be either a vector of classes (numeric, character or factor) or a two-column matrix with the counts for each class.")
  }

  # Handle NAs in labels
  if(info_type == 'long'){
    ix_na <- is.na(labels)
  } else if(info_type == 'summarized'){
    ix_na <- apply(labels, 1, function(x) any(is.na(x)))
  }
  if(any(ix_na)){
    if(all(ix_na)){
      stop('All the labels are missing (NA).')
    }
    predictions <- predictions[!ix_na]
    warning(sprintf('%d labels are missing (NA) and will be ignored. Proceeding with %d observations.', sum(ix_na), sum(1-ix_na)))
  }
  if(info_type == 'long'){
    labels <- labels[!ix_na]
  } else if(info_type == 'summarized'){
    labels <- labels[!ix_na,]
  }

  # Handle NAs in predictions
  ix_na <- is.na(predictions)
  if(all(ix_na)){
    special.values <- -9999999
    predictions <- rep(special.values[1], length(predictions))
    warning(sprintf('All the predictions are missing (NA) and will be replaced with the special value %d.', special.values[1]))
  } else if(any(ix_na)){
    special.values <- c(-9999999, special.values)
    predictions[ix_na] <- special.values[1]
    warning(sprintf('%d predictions are missing (NA) and will be replaced with the special value %d.', sum(ix_na), special.values[1]))
  }


  # Initial info
  out <- list()
  out$call <- match.call()
  class(out) <- 'wroc'

  # Homogeneous summarized format
  if(info_type == 'long'){
    predictions <- predictions
    labels <- data.frame(
      n_bad = as.numeric(labels == level.bad),
      n_good = as.numeric(labels != level.bad)
    )
  } else{
    predictions <- predictions
    labels <- as.matrix(labels)
    labels <- data.frame(
      n_bad = labels[,col.bad],
      n_good = labels[,which(!(1:2 %in% col.bad))]
    )
  }

  # Deal with special values separately
  if(!is.null(special.values)){
    special.values <- unique(special.values)
    special <- list()
    special_ix <- sapply(predictions, function(i){
      (i %in% special.values)
    })
    special$buckets <- as.numeric(sapply(predictions[special_ix], function(i){
      -which(i == special.values)
    }))
    special$predictions <- predictions[special_ix]
    special$n_bad <- labels[special_ix,'n_bad']
    special$n_good <- labels[special_ix,'n_good']
    special <- special %>%
      as.data.frame() %>%
      dplyr::group_by(buckets, predictions) %>%
      dplyr::summarise(n_bad = sum(n_bad),
                n_good = sum(n_good))

    regular <- list()
    regular$predictions <- predictions[!special_ix]
    regular$n_bad <- labels[!special_ix,'n_bad']
    regular$n_good <- labels[!special_ix,'n_good']
  } else {
    special_ix <- logical(0)
    special <- list()
    special$buckets <- NULL
    special$predictions <- NULL
    special$n_bad <- NULL
    special$n_good <- NULL

    regular <- list()
    regular$predictions <- predictions
    regular$n_bad <- labels$n_bad
    regular$n_good <- labels$n_good
  }


  # Initial grouping using qcut
  if(is.null(ngroups)){
    message('Using exact ROC curve. This may be very slow for continuous variables! Try using a smaller number for ngroups.')
    regular$buckets <- factor(regular$predictions)
  } else{
    # regular$buckets <- as.numeric(cut2(regular$predictions, g = ngroups))
    regular$buckets <- as.numeric(qcut3(regular$predictions,
                                        wt = regular$n_bad + regular$n_good,
                                        g = ngroups)$x)
  }

  # And summarizing
  ix <- c((special$n_bad + special$n_good != 0),
          (regular$n_bad + regular$n_good != 0))
  info <- data.frame(bucket = c(special$buckets, regular$buckets),
                     x = c(special$predictions, regular$predictions),
                     n_bad = c(special$n_bad, regular$n_bad),
                     n_good = c(special$n_good, regular$n_good)) %>%
    dplyr::filter(ix) %>%
    dplyr::group_by(bucket) %>%
    dplyr::summarise(population = sum(n_bad) + sum(n_good),
              lower_limit = max(x),
              upper_limit = max(x), # Gets lagged later
              n_bad = sum(n_bad))

  # Special population index
  ix <- (info$bucket %in% unique(special$buckets))

  # Totals
  totals <- list()
  totals$population <- sum(info$population)
  totals$spec_population <- sum(info$population[ix])
  totals$bad <- sum(info$n_bad)
  totals$spec_bad <- sum(info$n_bad[ix])
  totals$good <- totals$population - totals$bad
  totals$spec_good <- totals$spec_population - totals$spec_bad

  ## Finish info table
  out$info <- info %>%
    dplyr::mutate(n_good = population - n_bad,
           p_bad = n_bad / population,
           p_good = 1 - p_bad,
           ac_population = c(rep(NA, sum(ix)), cumsum(population[!ix])),
           ac_bad = c(rep(NA, sum(ix)), cumsum(n_bad[!ix])),
           ac_good = c(rep(NA, sum(ix)), cumsum(n_good[!ix])),
           d_population = population/totals$population)

  # Special cases if there are no obs of a class
  if(totals$bad != 0){
    aux <- out$info$n_bad/totals$bad
    if(totals$good != 0){
      out$info$d_bad <- aux
      out$info$d_good <- out$info$n_good/totals$good
    } else{
      out$info$d_bad <- aux
      out$info$d_good <- aux
    }
  } else if(totals$good != 0){
    aux <- out$info$n_good/totals$good
    out$info$d_bad <- aux
    out$info$d_good <- aux
  } else{
    stop('There are no goods or bads. Check your data.')
  }

  # Special cases if there are only special values for goods, bads or population
  if(totals$population != totals$spec_population){
    out$info$d_ac_population <- out$info$ac_population/(totals$population - totals$spec_population)

    if(totals$bad != totals$spec_bad){
      aux <- out$info$ac_bad/(totals$bad - totals$spec_bad)
      if(totals$good != totals$spec_good){
        out$info$d_ac_bad <- aux
        out$info$d_ac_good <- out$info$ac_good/(totals$good - totals$spec_good)
      } else{
        out$info$d_ac_bad <- aux
        out$info$d_ac_good <- aux
      }
    } else{
      if(totals$good != totals$spec_good){ # For clarity
        aux <- out$info$ac_good/(totals$good - totals$spec_good)
        out$info$d_ac_bad <- aux
        out$info$d_ac_good <- aux
      } else{
        stop("This can't happen.")
      }
    }

  } else{
    aux <- ifelse(ix, NA, 1)
    out$info$d_ac_population <- aux
    out$info$d_ac_bad <- aux
    out$info$d_ac_good <- aux
  }

  # Weight of Evidence
  out$info$woe <- ifelse(out$info$p_bad == 0 | out$info$p_good == 0,
                         NA,
                         log(out$info$p_good/out$info$p_bad))

  if(any(is.na(out$info$woe))){
    if(all(is.na(out$info$woe))){
      warning('All buckets have observations of a single class. Replacing WoE with +5 if there are only goods and -5 if there are only bads.')
      out$info <- out$info %>%
        dplyr::mutate(woe = ifelse(p_bad == 0, 5, -5))
    } else{
      warning('Some buckets have observations of a single class. Replacing WoE with twice the maximum / minimum WoE among other buckets.')
    out$info <- out$info %>%
      dplyr::mutate(woe = ifelse(is.na(woe),
                          ifelse(p_bad == 0,
                                 2*max(woe, na.rm = T),
                                 2*min(woe, na.rm = T)),
                          woe))
    }
  }

  out$special <- out$info[ix,] %>% intcols2double # Homogenize numeric types to double (for comparison purposes)

  out$info <- out$info[!ix,] %>%
    rbind(0, .)
  out$info$upper_limit[1] <- -Inf
  out$info$upper_limit[nrow(out$info)] <- Inf
  out$info$lower_limit <- lag(out$info$upper_limit)
  out$info$lower_limit[1] <- -Inf


  out$ngroups <- nrow(out$info) - 1
  out$nspecial <- nrow(out$special)
  out$totals <- totals

  if(length(special_ix > 0) && all(special_ix)){
    # This is a dummy
    suppressWarnings(
      suppressMessages(
        out$info <- wroc.default(1:3, c(0,1,1))$info[c(1,1),]
      )
    )
    out$info[2,c('bucket','upper_limit','d_ac_population','d_ac_bad','d_ac_good')] <- c(1, Inf, 1, 1, 1)
  }

  out
}

#' @describeIn wroc Formula-based for simple processing of several variables at
#'   once.
#' @export
wroc.formula <- function(formula, data, ngroups = NULL, level.bad=1,
                         special.values=NULL, verbose = (nrow(data) > 10000)){
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

  cat(sprintf('~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~\nGenerating %d ROC curves\nPrecision parameter (ngroups): %d\n\n',
              ncol(ds), ngroups))
  out <- lapply(1:ncol(ds), function(j){
    if(verbose){
      cat(sprintf('(%d/%d) Generating ROC curve for variable\t%s\t@ %s\n',
                  j, ncol(ds), names(ds)[j], Sys.time()))
    }
    if(reuse_special){
      spvals <- special.values
    } else {
      spvals <- special.values[[names(ds)[j]]]
    }
    wroc.default(ds[[j]], y, ngroups = ngroups, level.bad = level.bad,
         special.values = spvals)
  })
  cat(sprintf('\nFinished generating ROC curves @ %s\n',
              Sys.time()))
  names(out) <- names(ds)
  class(out) <- c('wroc.list','list')
  out
}

#' @rdname wroc
#' @export
print.wroc <- function(x){
  print(summary(x, performance = F), extended = T)
}

#' @rdname wroc
#' @export
c.wroc <- function(...){
  ar <- list(...)
  if(is.null(names(ar))){
    stop('The arguments to c.wroc must be named and the names must coincide with the names of the variables used to construct the wroc objects.')
  }
  class(ar) <- c('wroc.list', 'list')
  ar
}

#' @rdname wroc
#' @export
all.equal.wroc <- function(target, current, identical=FALSE){
  if(!identical){
    black.list <- 'call'
    if(nrow(target$special) == 0 && nrow(current$special) == 0){
      black.list <- c(black.list, 'special')
    }
    if(nrow(target$info) == 0 && nrow(current$info) == 0){
      black.list <- c(black.list, 'info')
    }
    target <- target[!(names(target) %in% black.list)]
    current <- current[!(names(current) %in% black.list)]
  }
  all.equal.list(target, current)
}


#' Plot ROC Curves
#'
#' Uses \code{ggplot2} to plot several interesting aspects about a ROC curve.
#'
#' @param x An object of class \code{wroc} or \code{wroc.list}.
#' @param type One of 'accum' (ROC using fnr ~ tnr), 'roc' (standard ROC tpr ~
#'   fpr), 'default' (probability of being in the positive class) and 'woe'
#'   (Weight of Evidence: log-odds of being in the *negative* class).
#' @param include.special Should special values be included in the plot? Only
#'   applies to 'default' and 'woe' options
#' @return An object of class ggplot which can be then be modified if needed. If
#'   \code{x} is a \code{wroc.list}, then a list of \code{ggplot} objects is
#'   returned. Depending on the value of \code{save.pdf}, a PDF with the plots
#'   is saved to \code{file} or the plots are shown one by one.
#' @export
plot.wroc <- function(x,
                      type = c('accum','roc','default','woe'),
                      include.special = TRUE,
                      label.size = 3){
  require(ggplot2)

  if(type[1] == 'accum'){
    p <- x$info %>%
      ggplot(aes(d_ac_good, d_ac_bad, color=bucket)) +
      geom_rect(xmin=0,xmax=1,ymin=0,ymax=1,fill=NA, color='black') +
      geom_segment(aes(x=0,xend=1,y=0,yend=1), color='black', linetype='dashed') +
      geom_line(size=1) +
      geom_point() +
      geom_point(size=1.5, shape=1, color='black') +
      scale_color_gradientn(colours = c('blue','green','yellow','red')) +
      labs(title='Accumulation', x='Cumulative % Goods (tnr)',y='Cumulative % Bads (fnr)') +
      coord_equal()
  } else if(type[1] == 'roc'){
    p <- x$info %>%
      ggplot(aes(1-d_ac_good, 1-d_ac_bad, color=bucket)) +
      geom_rect(xmin=0,xmax=1,ymin=0,ymax=1,fill=NA, color='black') +
      geom_segment(aes(x=0,xend=1,y=0,yend=1), color='black', linetype='dashed') +
      geom_line(size=1) +
      geom_point() +
      geom_point(size=1.5, shape=1, color='black') +
      scale_color_gradientn(colours = c('blue','green','yellow','red')) +
      labs(title='ROC Curve', x='fpr',y='tpr') +
      coord_equal()
  } else if(type[1] == 'default'){
    ds <- x$info[-1,]
    if(include.special) ds <- rbind(x$special, ds)
    ds <- ds %>%
      dplyr::mutate(barcol = ifelse(bucket < 0, 'salmon', 'darkgrey'),
             pointcol = ifelse(bucket < 0, 'salmon', 'black'))
    brks <- 1:nrow(ds)
    labls <- sprintf('B%d: (%.2f, %.2f]',
                     ds$bucket,
                     ds$lower_limit,
                     ds$upper_limit)
    if(include.special) labls[1:nrow(x$special)] <- gsub('\\(','[',labls[1:nrow(x$special)])
    labls[length(labls)] <- gsub(']',')',labls[length(labls)])
    p <- ds %>%
      dplyr::mutate(i = row_number(),
             norm_p_bad = p_bad*max(d_population)/max(p_bad)) %>%
      ggplot(aes(i)) +
      geom_bar(aes(y=d_population, fill=barcol), stat='identity') +
      geom_line(aes(y=norm_p_bad, color='black', group=pointcol)) +
      geom_point(aes(y=norm_p_bad, color=pointcol)) +
      # geom_point(y=norm_p_bad, color='black', shape=1) +
      geom_text(aes(y = norm_p_bad, label=sprintf('%.2f %%',100*norm_p_bad)),
                size = label.size, vjust=-1) +
      scale_color_identity() +
      scale_fill_identity() +
      scale_x_continuous(breaks=brks,labels=labls) +
      scale_y_continuous(labels = scales::percent) +
      theme(axis.text.x = element_text(angle=90)) +
      labs(title = 'Default Rate',
           x = 'Bucket',
           y = '% Population')

  } else if(type[1] == 'woe'){
    ds <- x$info[-1,]
    if(include.special) ds <- rbind(x$special, ds)
    ds <- ds %>%
      dplyr::mutate(barcol = ifelse(bucket < 0, 'salmon', 'darkgrey'),
             pointcol = ifelse(bucket < 0, 'salmon', 'black'))
    brks <- 1:nrow(ds)
    labls <- sprintf('B%d: (%.2f, %.2f]',
                     ds$bucket,
                     ds$lower_limit,
                     ds$upper_limit)
    if(include.special) labls[1:nrow(x$special)] <- gsub('\\(','[',labls[1:nrow(x$special)])
    labls[length(labls)] <- gsub(']',')',labls[length(labls)])
    p <- ds %>%
      dplyr::mutate(i = row_number(),
                    norm_woe = woe*max(d_population)/max(woe)) %>%
      ggplot(aes(i)) +
      geom_bar(aes(y=d_population, fill=barcol), stat='identity') +
      geom_line(aes(y=norm_woe, color='black', group=pointcol)) +
      geom_point(aes(y=norm_woe, color=pointcol)) +
      # geom_point(y=norm_p_bad, color='black', shape=1) +
      geom_text(aes(y = norm_woe, label=round(woe, 2)),
                size = label.size, vjust=-1) +
      scale_color_identity() +
      scale_fill_identity() +
      scale_x_continuous(breaks=brks,labels=labls) +
      scale_y_continuous(labels = scales::percent) +
      theme(axis.text.x = element_text(angle=90)) +
      labs(title = 'Weight of Evidence',
           x = 'Bucket',
           y = '% Population')
  }
  p
}
###############################################################


#' Optimal Bucketing With a Fixed Number of Buckets
#'
#' Optimizes the bucketing to have the maximum Gini index possible for the given
#' trend and number of buckets.
#'
#' @param x An object of class \code{wroc}.
#' @param nbuckets The number of buckets the final bucketing should have.
#' @param trend The trend according to the cummulative ROC curve.
#' @return An object of class \code{wroc} but with the original bucket numbers.
#'   The class \code{analyze.wroc} was added to show the curve has been
#'   modified. If the input is a \code{wroc.list} then the output is again a
#'   \code{wroc.list} of constraint-optimized ROC curves.
#' @export
analyze <- function(x, ...) UseMethod("analyze")

#' @describeIn analyze Method for a single variable.
#' @export
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
  out$gini_at_start <- gini_at_start
  out$gini <- gini
  out$gini_loss <- gini_loss
  out$class <- c('analyze.wroc', class(out))
  out
}

###############################################################

#' Optimal Gini Bucketing
#'
#' Optimizes the bucketing to have the maximum Gini index possible for the given
#' trend.
#'
#' @param x An object of class \code{wroc}.
#'
#' @return An object of class \code{wroc}. The class \code{optimize.wroc} has
#'   been added to show the fact that the original curve has been optimized. If
#'   the input is a \code{wroc.list} then the output is again a \code{wroc.list}
#'   of optimized ROC curves.
#' @seealso \code{analyze.wroc}
#' @export
optimize <- function(x, ...) UseMethod("optimize")

# Internal function for optimize.wroc (optimal algorithm)
optimize_optimal_ <- function(x, trend = c('auto','upper','lower')){
  ds <- x$info
  if(trend[1] == 'auto'){
    trend <- choose.trend_(x)
  }

  fun <- ifelse(trend[1] == 'upper', max, min)

  keepers <- c(1)
  i <- 1
  while(i <= nrow(ds)-1){
    dx <- ds$d_ac_good[(i+1):nrow(ds)] - ds$d_ac_good[i]
    dy <- ds$d_ac_bad[(i+1):nrow(ds)] - ds$d_ac_bad[i]
    slopes <- ifelse(dx == 0, Inf, dy/dx)
    extrema_ix <- which(slopes == fun(slopes))
    j <- max(extrema_ix)
    keepers <- c(keepers, i + j)
    i <- i + j
  }

  buckets_to_remove <- ds$bucket[!(1:nrow(ds) %in% keepers)]
  out <- subset(x, buckets = buckets_to_remove)

  # out$removed.buckets <- c(x$removed.buckets, buckets_to_remove)
  out$trend <- trend
  out$call.optimize <- match.call()

  class(out) <- c('optimal.wroc', 'wroc')
  out
}

# Internal function for optimize.wroc (optimal preserving population)
optimize_magic_ <- function(x, trend = c('auto','upper','lower'), min_p_pob=0.5){
  ds <- x$info
  if(trend[1] == 'auto'){
    trend <- choose.trend_(x)
  }
  if(min_p_pob > 1 || min_p_pob < 0){
    stop('min_p_pob must be between 0 and 1.')
  }

  # boolean for 'decreasing' flag in order() call
  decreasing <- (trend[1] == 'upper')

  keepers <- c(1)
  i <- 1
  while(i <= nrow(ds)-1){
    dx <- ds$d_ac_good[(i+1):nrow(ds)] - ds$d_ac_good[i]
    dy <- ds$d_ac_bad[(i+1):nrow(ds)] - ds$d_ac_bad[i]
    slopes <- ifelse(dx == 0, Inf, dy/dx)
    or <- order(slopes, decreasing = decreasing)
    len <- length(or)
    for(k in 1:len){
      j <- or[k]
      if(ds$d_ac_population[i+j] - ds$d_ac_population[i] >= min_p_pob ||
         k == len){
        if(ds$d_ac_population[i+j] > (1 - min_p_pob)){
          j <- len
        }
        keepers <- c(keepers, i + j)
        i <- i + j
        break
      }
    }
  }

  buckets_to_remove <- ds$bucket[!(1:nrow(ds) %in% keepers)]
  out <- subset(x, buckets = buckets_to_remove)

  # out$removed.buckets <- c(x$removed.buckets, buckets_to_remove)
  out$trend <- trend
  out$call.optimize <- match.call()

  class(out) <- c('magic.wroc', 'wroc')
  out
}

#' @describeIn optimize Optimize a single ROC curve.
#' @inheritParams analyze.wroc
#' @param method One of "optimal" (optimal Gini regardless of population),
#'   "magic" [default] (tries to optimize Gini while keeping the distribution as
#'   homogeneous as possible) and "cascade" (iteratively pastes the smallest
#'   bucket left or right so as to lose as little Gini as possible).
#' @param min_p_pob Numeric between 0 and 1. Minimum proportion of population to
#'   keep in each bucket. Only for methods "magic" and "cascade",
#' @export
optimize.wroc <- function(x, trend = c('auto','upper','lower'), method = c('magic', 'optimal', 'cascade'), min_p_pob = 0.05){
  if('optimal.wroc' %in% class(x)){
    warning('This is an already optimal ROC curve. Returning the input ROC curve. ')
    return(x)
  } else if('subset.wroc' %in% class(x)){
    warning('This ROC curve has been manually tampered with. Resetting bucket names.')
    x <- reset.buckets(x)
  }
  if(!(trend[1] %in% c('auto','upper','lower'))){
    warning('Invalid trend. Only "upper", "lower" or "auto" are valid. Defaulting to "auto".')
    trend <- 'auto'
  }
  if(trend[1] == 'auto'){
    trend <- choose.trend_(x)
  }
  if(method[1] == 'optimal'){
    optimize_optimal_(x, trend)
  } else if(method[1] == 'magic'){
    optimize_magic_(x, trend, min_p_pob = min_p_pob)
  } else if(method[1] == 'cascade'){
    warning('Cascade method is not yet available.')
  } else{
    stop('Invalid method. Please choose one of "optimal", "magic" or "cascade".')
  }
}

###############################################################

#' Performance Metrics for ROC Curves
#'
#' Calculates several performance metrics for ROC curves, including AUC and Gini
#' index for both trends.
#'
#' @param x An object of type \code{wroc} or \code{wroc.list}.
#' @return A list with several performance statistics *based on the inverse ROC
#'   curve* (as given by the 'trend' option in \code{plot.wroc}), that is fnr ~
#'   tnr:
#'
#'   \itemize{ \item{\code{auc}: Area under the ROC curve}
#'   \item{\code{gini_two_way}: abs(2*auc - 1) (different from normal Gini when
#'   the ROC curves crosses the identity line)} \item{\code{gini_up}: Gini with
#'   an upper trend} \item{\code{gini_down}: Gini with a lower trend}
#'   \item{\code{gini}: The larger of gini_up and gini_down}
#'   \item{\code{best_trend}: The trend of the larger Gini} \item{\code{ks}:
#'   Kolmogorov-Smirnov statistic} \item{\code{iv}: Information Value} } In the
#'   case where \code{x} is a \code{wroc.list}, then a list of the
#'   aforementioned is returned instead.
#' @export
performance <- function(x, ...) UseMethod("performance")

#' @describeIn performance Performance measures of a \code{wroc} object.
#' @export
performance.wroc <- function(x){
  ds <- x$info

  if(x$ngroups == 0){ # There are only special values
    tr <- data.frame(
      area_auc = 0.5,
      area_gini_up = 0,
      area_gini_down = 0,
      ks_prospect = 0,
      iv_contribution = 0
    )
  } else{
    tr <- ds %>%
      dplyr::mutate(dx = (d_ac_good - dplyr::lag(d_ac_good)),
             bases_auc = (d_ac_bad + dplyr::lag(d_ac_bad)),
             bases_gini_up_raw = pmax(d_ac_bad - d_ac_good, 0),
             bases_gini_up = bases_gini_up_raw + lag(bases_gini_up_raw),
             bases_gini_down_raw = pmax(d_ac_good - d_ac_bad, 0),
             bases_gini_down = bases_gini_down_raw + lag(bases_gini_down_raw),
             area_auc = bases_auc*dx/2,
             area_gini_up = bases_gini_up*dx/2,
             area_gini_down = bases_gini_down*dx/2,
             ks_prospect = abs(d_ac_good - d_ac_bad),
             iv_contribution = (d_good - d_bad)*woe) %>%
      .[-1,]
  }

  out <- list()
  out$auc <- sum(tr$area_auc)
  out$gini_two_way <- abs(2*out$auc - 1)
  out$gini_up <- 2*sum(tr$area_gini_up)
  out$gini_down <- 2*sum(tr$area_gini_down)
  out$gini <- max(out$gini_up, out$gini_down)
  out$best_trend <- ifelse(out$gini_up > out$gini_down, 'upper', 'lower')
  out$ks <- max(tr$ks_prospect)
  out$iv <- sum(tr$iv_contribution)

  class(out) <- c('wroc.performance', 'list')
  out
}

#' @describeIn performance Quick evaluation without need of a \code{wroc}
#'   object.
#' @inheritParams wroc
#' @param ... Aditional arguments for \code{wroc}
#' @export
performance.numeric <- function(predictions, labels, ...){
  suppressWarnings(x <- wroc.default(predictions, labels, ...))
  performance.wroc(x)
}

###############################################################

#' Transform Variables with a ROC Curve
#'
#' Transforms raw variables into bucketed ones. The pasted value can be any
#' column in the \code{object$info} \code{data.frame}, but it usually is the
#' WoE, the bucket number or the probability of being in the positive class.
#'
#' @param object An object of class \code{wroc} or \code{wroc.list}.
#' @param newdata A numeric vector containing values of the raw variable to be
#'   transformed. If \code{object} is a \code{wroc.list}, then \code{newdata}
#'   must be a \code{data.frame} with the variables in the list.
#' @param type A string containing the name of the column in \code{object$info},
#'   to be used as transformation, usually 'woe', 'bucket' or 'p_bad'.
#' @return A vector of the same length as \code{newdata} containing the
#'   transformed variable. If \code{object} is a \code{wroc.list}, then the
#'   output is a \code{data.frame} of the transformed variables. In that case,
#'   \code{keep.data} can be used to append the transformed variables to the
#'   original dataset instead of just returning them.
#' @export
predict.wroc  <- function(object,
                          newdata,
                          type = c('woe','bucket','p_bad')){

  type <- type[1]
  if(!(type %in% names(object$info))){
    stop('Type must be a column name of object$info, such as "woe", "p_bad" and "bucket".')
  }
  if(!is.numeric(newdata)){
    warning('newdata must be numeric. Attempting to coerce to numeric.')
    newdata <- as.numeric(newdata)
  }

  out <- rep(NA, length(newdata))
  special_ix <- (newdata %in% object$special$lower_limit)

  out[!special_ix] <- cut(x = newdata[!special_ix],
                          breaks = c(-Inf, object$info$upper_limit[-1]),
                          include.lowest = FALSE,
                          labels = object$info$bucket[-1],
                          right = TRUE) %>%
    as.character() %>%
    as.numeric()
  out[special_ix] <- newdata[special_ix]

  vals <- rbind(
    data.frame(
      id = object$info$bucket[-1],
      type = object$info[[type]][-1]
    ),
    data.frame(
      id = object$special$upper_limit,
      type = object$special[[type]]
    )
  )
  yhat <- data.frame(id=out) %>%
    left_join(vals, by='id') %>%
    .$type
  yhat
}

###############################################################

#' Paste Adjacent Buckets of a ROC Curve
#'
#' Pastes the given bucket names to the immediately following one. This allows
#' for manual tampering of the ROC curves.
#'
#' @param x An object of class \code{wroc}.
#' @param buckets Bucket numbers to be pasted to the right.
#' @return An object of class \code{wroc} with the given buckets pasted to the
#'   right. The class \code{subset.wroc} has been added to show that the ROC
#'   curve isn't the original one. The original bucket numbers are used for
#'   clarity (you can use reset.buckets afterwards to change this)
#' @export
subset.wroc <- function(x, buckets = NULL, ...){
  if(is.null(buckets) || length(buckets) == 0) return(x)
  if(any(buckets < 0)){
    stop('Cannot paste special valued buckets. Join the levels by hand before calling wroc.')
  } else {
    ds <- x$info[-1,]
    sp <- x$special
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

  suppressWarnings(
    out <- wroc(predictions=c(sp$bucket, ds$bucket),
                labels=rbind(
                  cbind(sp$n_bad, sp$n_good),
                  cbind(ds$n_bad, ds$n_good)
                ),
                ngroups = NULL,
                special.values = x$special$bucket)
  )


  out$special <- x$special
  out$nspecial <- x$nspecial
  out$info$bucket <- c(0, ds$bucket[-ix])
  out$info$lower_limit <- c(-Inf, ds$lower_limit[-(ix+1)])
  out$info$upper_limit <- c(-Inf, ds$upper_limit[-(ix)])
  ### out$ngroups <- nrow(out$info) - 1
  out$call.subset <- match.call()
  out$pasted.buckets <- buckets
  class(out) <- c('subset.wroc', class(out))
  out
}

###############################################################

#' Important Information About ROC Curves
#'
#' Returns and prints basic important information about a \code{wroc} or
#' \code{wroc.list} object.
#'
#' @param object An object of class \code{wroc}.
#' @param performance Calculate performance measures as given by
#'   \code{performance} (AUC, Gini, KS, IV)
#' @return An object of class \code{summary.wroc}, which has interesting
#'   non-internal information about the ROC curve and prints prettily. If
#'   performance = FALSE then performance is not run and hence the method will
#'   be faster for ROC curves with a very large number of buckets. If
#'   \code{object} is a list, then the return object has both a list of the
#'   individual summaries and a table showing the most important information
#'   about them in a compact way.
#' @export
summary.wroc <- function(object, performance = TRUE, ...){
  out <- list()
  if(object$nspecial > 0){
    spvals <- cbind(data.frame(type='special', stringsAsFactors = F),
                    object$special)
  } else{
    spvals <- NULL
  }

  out$info <- rbind(
      spvals,
      cbind(data.frame(type='normal', stringsAsFactors = F), object$info[-1,])
    ) %>%
    dplyr::mutate(range = ifelse(row_number() < n(),
                          sprintf('(%.2f, %.2f]', lower_limit, upper_limit),
                          sprintf('(%.2f, %.2f)', lower_limit, upper_limit))) %>%
    dplyr::select(bucket, type, lower_limit, upper_limit, range,
                  n_good, n_bad, population, d_good, d_bad, d_population,
                  d_ac_good, d_ac_bad, d_ac_population, p_bad, woe) %>%
    dplyr::as.tbl()
  if(performance){
    out <- c(out, performance.wroc(object))
  }
  out$performance <- performance
  out$totals <- object$totals
  class(out) <- c('summary.wroc')
  out
}

#' @rdname summary.wroc
#' @param extended Should details be printed?
#' @export
print.summary.wroc <- function(object, extended = FALSE, ...){
  regular <- dplyr::filter(object$info, type == 'normal')
  special <- dplyr::filter(object$info, type == 'special')
  tot <- object$totals
  out <- sprintf(
    'Regular Buckets: %d\nSpecial Buckets: %d\nTotal Population: %d\nSpecial Population: %d (%.2f%%)\nSmallest bucket Pop.: %d (%.2f%%)',
    nrow(regular), nrow(special),
    tot$population,
    tot$spec_population,
    100*tot$spec_population/tot$population,
    min(regular$population),
    100*min(regular$d_population)
  )

  if(object$performance){
    out <- sprintf(
      '%s\n\nGini Index: %.2f%% (trend: %s)\nAUC: %.3f\nKolmogorov-Smirnov Statistic: %.3f\nInformation Value: %.3f)',
      out,
      100*object$gini, object$best_trend, object$auc, object$ks, object$iv
    )
  }
  cat(out)
  if(extended){
    cat('\n\nDetails:\n\n')
    print(object$info)
  }
}


###############################################################

#' Export WROCs and WROC Lists
#'
#' Serialize \code{wroc}'s, \code{wroc.list}'s or a list summaries as tables. There is
#' also an option for copying them to clipboard for easy pasting into Excel.
#'
#' @param x An object of class \code{wroc} or \code{wroc.list}.
#' @param export Should the table be exported? (As opposed to copied to the
#'   clipboard).
#' @param ... Additional parameters for \code{write.table}.
#' @export
copy <- function(x, ...) UseMethod("copy")

#' @rdname copy
#' @export
copy.wroc <- function(x, export=FALSE, ...){
  ellipsis <- list(...)
  if(export){
    if(is.null(ellipsis$file)){
      file <- sprintf('WROC List %s', Sys.time())
    }
  } else{
    if(is.null(ellipsis$file)){
      file <- 'clipboard-12345'
    }
  }

  temp <- summary(x)$info
  temp$lower_limit <- gsub('-Inf', 'minus Inf', as.character(temp$lower_limit))
  if(export){
    write.table(temp, row.names = F, file = file, ...)
  } else{
    write.table(temp, file = file, sep = '\t', row.names = F, ...)
  }
}

#' @rdname copy
#' @export
copy.wroc.list <- function(x, export=FALSE, ...){
  ellipsis <- list(...)
  if(export){
    if(is.null(ellipsis$file)){
      file <- sprintf('WROC List %s', Sys.time())
    }
  } else{
    if(is.null(ellipsis$file)){
      file <- 'clipboard-12345'
    }
  }

  temp <- lapply(1:length(x), function(i){
    cbind(data.frame(variable=names(x)[i], stringsAsFactors = F),
          summary(x[[i]])$info,
          stringsAsFactors = F)
  }) %>%
    rbind_all
  temp$lower_limit <- gsub('-Inf', 'minus Inf', as.character(temp$lower_limit))
  if(export){
    write.table(temp, row.names = F, file = file, ...)
  } else{
    write.table(temp, row.names = F, file = file, sep = '\t', ...)
  }
}

#' @rdname copy
#' @export
copy.summary.wroc.list <-  function(x, export=FALSE, ...){
  if(export){
    write.table(x$info, row.names = F, file = file, ...)
  } else{
    write.table(x$info, row.names = F, file='clipboard-12345', sep = '\t', ...)
  }
}



#compactify.wroc






