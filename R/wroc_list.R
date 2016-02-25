
require(dplyr)
require(tidyr)
require(ggplot2)

# wroc.list methods
`[.wroc.list` <- function(x, i){
  out <- x
  class(out) <- 'list'
  out <- out[i]
  class(out) <- c('wroc.list','list')
  out
}

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

optimize.wroc.list <- function(x, trends = 'auto'){
  if((length(trends) != length(x)) && (trends != 'auto')){
    stop('trends must be either "auto" or a character vector with a trend for each vriable.')
  } else if((length(trends) == 1) && (trends == 'auto')){
    trends <- rep('auto', length(x))
  }
  for(i in 1:length(x)){
    x[[i]] <- optimize.wroc(x[[i]], trends[i])
  }
  x
}

predict.wroc.list <- function(object,
                              newdata,
                              type = c('woe','bucket','p_bad'),
                              keep.data = FALSE,
                              prefix = type){
  yhats <- list()
  for(i in 1:length(object)){
    variable <- names(object)[i]
    yhats[[variable]] <- predict(object[[i]],
                                 newdata,
                                 variable,
                                 type,
                                 keep.data = FALSE,
                                 prefix)
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
