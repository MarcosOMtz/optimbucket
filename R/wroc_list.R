
require(dplyr)
require(tidyr)
require(ggplot2)

#' @export
print.wroc.list <- function(x, ...){
  for(i in 1:length(x)){
    if(i > 1) cat(sprintf('\n\n'))
    cat(sprintf('=================================\nVariable: %s\n\n', names(x)[i]))
    print(x[[i]])
  }
}

#' @export
`[.wroc.list` <- function(x, i){
  out <- x
  class(out) <- 'list'
  out <- out[i]
  class(out) <- c('wroc.list','list')
  out
}

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

#' @describeIn optimize Runs the algorithm on a \code{wroc.list}.
#' @inheritParams analyze.wroc
#' @param trends A character vector containing the trend to choose for each
#'   variable or 'auto' to let the algorithm find out on its own.
#' @export
optimize.wroc.list <- function(x, trends = 'auto', verbose = TRUE){
  if((length(trends) != length(x)) && (trends != 'auto')){
    stop('trends must be either "auto" or a character vector with a trend for each variable.')
  } else if((length(trends) == 1) && (trends == 'auto')){
    trends <- rep('auto', length(x))
  }
  for(i in 1:length(x)){
    if(verbose){
      cat(sprintf('Variable # %d (%.0f %%):\t%s\n', i, 100*i/length(x), names(x)[i]))
    }
    x[[i]] <- optimize.wroc(x[[i]], trends[i])
  }
  x
}

predict.wroc.list <- function(object,
                              newdata,
                              type = c('woe','bucket','p_bad'),
                              keep.data = FALSE,
                              prefix = type,
                              verbose = (nrow(newdata)*length(object) > 10000)){
  yhats <- list()
  for(i in 1:length(object)){
    if(verbose){
      cat(sprintf('Variable # %d (%.0f %%):\t%s\n',
                  i, 100*i/length(object), names(object)[i]))
    }
    variable <- names(object)[i]
    yhats[[variable]] <- predict(object[[i]],
                                 newdata[[variable]],
                                 type)
  }
  names(yhats) <- sprintf('%s_%s', prefix, names(object))
  yhats <- as.data.frame(yhats)
  if(keep.data){
    return(cbind(newdata, yhats))
  } else{
    return(yhats)
  }
}

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

summary.wroc.list <- function(x, performance = T, ...){
  out <- lapply(x, function(y){
    summary(y, performance)
  })
  class(out) <- c('summary.wroc.list', 'list')
  out
}

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

