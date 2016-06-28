
###############################################################

#' @rdname wroc
#' @export
print.wroc.list <- function(x, ...){
  for(i in 1:length(x)){
    if(i > 1) cat(sprintf('\n\n'))
    cat(sprintf('=================================\nVariable: %s\n\n', names(x)[i]))
    print(x[[i]])
  }
}

#' @rdname wroc
#' @export
c.wroc.list <- function(...){
  ar <- list(...)
  vars <- unlist(lapply(ar, names))
  if(any(duplicated(vars))){
    stop(sprintf('The following variables are present in two or more of the wroc.list objects: %s. Please make sure each variable appears only once.',
                 paste(unique(vars[duplicated(vars)]), collapse = ', ')))
  }
  out <- unlist(ar, recursive = FALSE)
  class(out) <- c('wroc.list', 'list')
  out
}

# #' @describeIn plot.wroc Plots for all variables in a list

#' @rdname plot.wroc
#' @param save.pdf Should the plots be saved to a PDF file or just shown?
#' @param file Name (with path) of the PDF file to output the plots to
#' @param ... Additional options to pass on to \code{pdf()}
#' @export
plot.wroc.list <- function(x,
                           type = c("accum", "roc", "default", "woe"),
                           include.special = TRUE,
                           save.pdf = FALSE,
                           file = 'Wroc Plots.pdf',
                           ...){
  if(save.pdf){
    pdf(file = file, ...)
  } else{
    cat(sprintf('~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~\nPlotting %s of a list of %d variables...', type[1], length(x)))
  }
  out <- list()
  for(i in 1:length(x)){
    if(type[1] == 'trend'){
      type <- 'default'
      message('type = "trend" will be deprecated in future versions. Use "default" instead')
    }
    if(type[1] == 'accum') ti <- sprintf('Accumulation: %s', names(x)[i])
    else if(type[1] == 'roc') ti <- sprintf('ROC: %s', names(x)[i])
    else if(type[1] == 'default') ti <- sprintf('Default Rate: %s', names(x)[i])
    else if(type[1] == 'woe') ti <- sprintf('WoE: %s', names(x)[i])
    else if(type[1] == 'points') ti <- sprintf('Points: %s', names(x)[i])
    print(
      out[[names(x)[i]]] <- plot(x[[i]],
                                 type = type,
                                 include.special = include.special) +
        labs(title=ti)
    )
    if(!save.pdf){
      readline(sprintf('(%d/%d) %s plot of variable %s',
                       i, length(x), type[1], names(x)[i]))
    }
  }
  if(save.pdf){
    dev.off()
  }
  out
}
###############################################################

#' @describeIn optimize Runs the algorithm on a \code{wroc.list}.
#' @inheritParams analyze
#' @param trends A character vector containing the trend to choose for each
#'   variable or 'auto' to let the algorithm find out on its own.
#' @export
optimize.wroc.list <- function(x, trends = 'auto', method = c('magic', 'optimal', 'cascade'), min_p_pob = 0.05, verbose = TRUE){
  if((length(trends) != length(x)) && (trends != 'auto')){
    stop('trends must be either "auto" or a character vector with a trend for each variable.')
  } else if((length(trends) == 1) && (trends == 'auto')){
    trends <- rep('auto', length(x))
  }
  cat(sprintf('~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~\nOptimizing %d variables...\n\n',
              length(x)))
  for(i in 1:length(x)){
    if(verbose){
      cat(sprintf('(%d/%d) Optimizing variable %s\t@ %s\n',
                  i, length(x), names(x)[i], Sys.time()))
    }
    x[[i]] <- optimize.wroc(x[[i]],
                            trend = trends[i],
                            method = method,
                            min_p_pob = min_p_pob)
  }
  cat(sprintf('\nFinished optimizing ROC curves @ %s\n',
              Sys.time()))
  x
}

###############################################################

#' @rdname predict.wroc
#' @inheritParams predict.wroc
#' @param keep.data Should predictions be appended to original data or should
#'   the predictions exclusively be returned?
#' @param prefix The prefix for the new variables. Defaults to \code{type}
#' @param verbose Should progress info be displayed?
#' @export
predict.wroc.list <- function(object,
                              newdata,
                              type = c('woe','bucket','p_bad'),
                              keep.data = FALSE,
                              prefix = type[1],
                              verbose = (nrow(newdata)*length(object) > 10000)){
  if(verbose){
    cat(sprintf('~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~\nTransforming %d variables to %s...\n\n',
              length(object), type[1]))
  }
  yhats <- list()
  for(i in 1:length(object)){
    if(verbose){
      cat(sprintf('(%d/%d) Transforming variable %s\t@ %s\n',
                  i, length(object), names(object)[i], Sys.time()))
    }
    variable <- names(object)[i]
    yhats[[variable]] <- predict(object[[i]],
                                 newdata[[variable]],
                                 type)
  }
  names(yhats) <- sprintf('%s_%s', prefix, names(object))
  yhats <- as.data.frame(yhats)
  cat(sprintf('\nFinished generating ROC curves @ %s\n',
              Sys.time()))
  if(keep.data){
    return(cbind(newdata, yhats))
  } else{
    return(yhats)
  }
}

###############################################################

#' Choose ROC Curves Manually from a List
#'
#' Provides a 'keep' and 'drop'-style for selecting ROC curves in a large list.
#'
#' @param x An object of type \code{wroc.list}
#' @param keep A character vector of the names of the ROC curves to be kept
#'   (overrides \code{drop})
#' @param drop A character vector of the names of the ROC curves to be discarded
#'   (is overriden by \code{keep})
#' @export
subset.wroc.list <- function(x, keep=NULL, drop=NULL){
  if((!(is.character(keep)) && !is.null(keep))
     ||
     (!(is.character(drop)) && !is.null(drop))){
    stop('keep and drop must be either character vectors or have a NULL value.')
  }
  if(is.null(keep)){
    keep <- names(x)
  } else {
    drop <- drop[which(!(drop %in% keep))]
  }

  ix <- which((names(x) %in% keep) & !(names(x) %in% drop))
  x[ix]
}

#' @rdname subset.wroc.list
#' @export
`[.wroc.list` <- function(x, i){
  out <- x
  class(out) <- 'list'
  out <- out[i]
  class(out) <- c('wroc.list','list')
  out
}

###############################################################

#' @describeIn performance Method for a list of ROC curves
#' @export
performance.wroc.list <- function(x, out.type = c('list', 'data.frame')){
  if(out.type[1] == 'data.frame'){
    out <- dplyr::rbind_all(lapply(1:length(x), function(i){
      data.frame(var=names(x)[i], as.list(performance(x[[i]])), stringsAsFactors = F)
    }))
  } else if(out.type[1] == 'list'){
    out <- lapply(x, function(y){
      as.list(performance(y))
    })
  } else{
    stop('out.type must be either "list" or "data.frame".')
  }
 return(out)
}

###############################################################

#' @rdname summary.wroc
#' @export
summary.wroc.list <- function(x, performance = T, ...){
  out <- list()
  out$summaries <- lapply(x, function(y){
    summary(y, performance)
  })
  variable <- names(out$summaries)
  tot_pop <- sapply(out$summaries, function(y) y$totals$population)
  spec_pop <- sapply(out$summaries, function(y) y$totals$spec_population)
  no_spec_pop <- tot_pop - spec_pop
  p_pop_no_spec <- 1 - (spec_pop/tot_pop)
  p_pop_spec <- 1 - p_pop_no_spec
  n_bad_no_spec <- sapply(out$summaries, function(y) y$totals$bad - y$totals$spec_bad)
  p_bad_no_spec <- n_bad_no_spec/(tot_pop - spec_pop)
  p_pop_smallest <- sapply(out$summaries, function(y){
    regular <- dplyr::filter(y$info, type == 'normal')
    min(regular$d_population)
  })
  gini <- sapply(out$summaries, function(y) y$gini)
  out$info <- data.frame(
    variable = variable,
    tot_pop = tot_pop,
    spec_pop = spec_pop,
    no_spec_pop = no_spec_pop,
    p_pop_no_spec = p_pop_no_spec,
    p_pop_spec = p_pop_spec,
    p_pop_smallest = p_pop_smallest,
    p_bad_no_spec = p_bad_no_spec,
    gini_no_spec = gini
  ) %>%
    dplyr::arrange(dplyr::desc(gini)) %>%
    cbind(rank=1:nrow(.), .)
  class(out) <- c('summary.wroc.list', 'list')
  out
}

#' @rdname summary.wroc
#' @export
print.summary.wroc.list <- function(x, ...){
  names_cond <- (mean(nchar(names(x$summaries))) > 15 | length(names(x$summaries)) > 10)
  if(names_cond){
    sep_names <- ',\n\t'
    tab_names <- '\t'
  } else{
    sep_names <- ', '
    tab_names <- ''
  }

  cat(paste(
    sprintf('~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~'),
    sprintf('List of %d ROC curves', length(x$summaries)),
    sprintf('Variables:'),
    paste0(tab_names, paste(names(x$summaries), collapse = sep_names)),
    sprintf('\nDetails:\n\n'),
    sep = '\n'))
  y <- x$info
  vars_num <- c('tot_pop','spec_pop','no_spec_pop')
  vars_perc <- grep('p_', names(y), value = T)
  y[vars_num] <- apply(y[vars_num], 2, function(s){
    format(s, scientific = F, big.mark = ',')
  })
  y[vars_perc] <- apply(y[vars_perc], 2, function(s){
    sprintf('%.2f%%', 100*s)
  })
  y$gini_no_spec <- round(y$gini_no_spec, 3)
  print(y)
}

