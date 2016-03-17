
require(dplyr)
require(tidyr)
require(ggplot2)

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
                           type = c("accum", "roc", "trend", "woe"),
                           include.special = TRUE,
                           save.pdf = TRUE,
                           file = 'Wroc Plots.pdf',
                           ...){
  if(save.pdf){
    pdf(file = file, ...)
  } else{
    cat(sprintf('~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~\nPlotting %s of a list of %d variables...', type[1], length(x)))
  }
  out <- list()
  for(i in 1:length(x)){
    print(
      out[[names(x)[i]]] <- plot(x[[i]],
                                 type = type,
                                 include.special = include.special) +
        labs(title=names(x)[i])
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
optimize.wroc.list <- function(x, trends = 'auto', verbose = TRUE){
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
    x[[i]] <- optimize.wroc(x[[i]], trends[i])
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
  cat(sprintf('~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~\nTransforming %d variables to %s...\n\n',
              length(object), type[1]))
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
    out <- rbind_all(lapply(1:length(x), function(i){
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
  out <- lapply(x, function(y){
    summary(y, performance)
  })
  class(out) <- c('summary.wroc.list', 'list')
  out
}

#' @rdname summary.wroc
#' @export
print.summary.wroc.list <- function(x, ...){
  names_cond <- (mean(nchar(names(x))) > 15 | length(names(x)) > 10)
  if(names_cond){
    sep_names <- ',\n\t'
    tab_names <- '\t'
  } else{
    sep_names <- ', '
    tab_names <- ''
  }

  cat(paste(
    sprintf('~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~'),
    sprintf('List of %d ROC curves', length(x)),
    sprintf('Variables:'),
    paste0(tab_names, paste(names(x), collapse = sep_names)),
    sprintf('\nDetails:\n\n'),
    sep = '\n'))
  for(i in 1:length(x)){
    if(i > 1) cat(sprintf('\n\n'))
    cat(sprintf('=================================\nVariable: %s\n\n', names(x)[i]))
    print(x[[i]])
  }
}

