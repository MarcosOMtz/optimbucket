
require(dplyr)
require(tidyr)
require(ggplot2)
require(testthat)
require(optimbucket)
set.seed(1234)


generate_sample_data <- function(ng=10000, nb=1000){

  goods <- data.frame(
    x = rnorm(ng, mean = ifelse(runif(ng) > 0.2, -2, 1)),
    y = 0
  )
  bads <- data.frame(
    x = rnorm(nb),
    y = 1
  )

  d <- rbind(goods, bads) %>%
    mutate(y_numeric = y,
           y_all_equal = 1,
           y_na = ifelse(runif(n()) < 0.5, NA, y),
           y_all_na = NA,
           y_factor = factor(y),
           y_character = as.character(y),
           y_multilevel = sample(letters, n(), replace=TRUE),
           x_factor = factor(x),
           x_all_equal = 3,
           x_factor2 = factor(round(10*x)),
           x_character = as.character(x),
           x_spec = ifelse(runif(n()) < 0.1,
                           -999,
                           ifelse(runif(n()) < 0.1,
                                  -998,
                                  x)),
           x_na = ifelse(runif(n()) < 0.5, NA, x),
           x_all_na = NA)
  d
}

# generate_sample_data(10,5)

# wroc

test_that('wroc can handle all sorts of labels', {
  d <- generate_sample_data(10000, 1000)

  expect_is(wroc(predictions=d$x, labels=d$y_numeric, ngroups = 50), 'wroc')
  expect_is(wroc(predictions=d$x, labels=d$y_factor, ngroups = 50), 'wroc')
  expect_is(wroc(predictions=d$x, labels=d$y_character, ngroups = 50), 'wroc')
  expect_is(wroc(predictions=d$x, labels=d$y_multilevel, ngroups = 50), 'wroc')
  expect_is(wroc(predictions=d$x, labels=d$y_all_equal, ngroups = 50), 'wroc')

  expect_warning(wroc(predictions=d$x, labels=d$y_na, ngroups = 50), 'NA|labels')
  expect_error(wroc(predictions=d$x, labels=d$y_all_na, ngroups = 50), 'NA|labels')
})

test_that('wroc can handle all sorts of predictions', {
  d <- generate_sample_data(10000, 1000)

  expect_is(wroc(predictions=d$x, labels=d$y, ngroups = 50), 'wroc')
  expect_is(wroc(predictions=d$x_factor, labels=d$y, ngroups = 50), 'wroc')
  expect_is(wroc(predictions=d$x_factor2, labels=d$y, ngroups = 50), 'wroc')
  expect_is(wroc(predictions=d$x_character, labels=d$y, ngroups = 50), 'wroc')

  expect_is(wroc(predictions=d$x_spec, labels=d$y, ngroups = 50, special.values = c(-999,-998)), 'wroc')
  expect_warning(wroc(predictions=d$x_na, labels=d$y, ngroups = 50), 'NA|predictions')
  expect_warning(wroc(predictions=d$x_all_na, labels=d$y, ngroups = 50), 'NA|predictions|All|all')
})

test_that("wroc's simple methods work",{
  d <- generate_sample_data(10000, 1000)
  w <- wroc(predictions=d$x_spec, labels=d$y, ngroups = 50,
            special.values = c(-999, -998))

  expect_null(print(w))
  expect_is(c(a=w, b=w), 'wroc.list')
  expect_warning(c(a=w, w), 'name')

  ### SUBSET
})
# Grouped mode
sum_data <- function(dat, x, y, sp=NULL, output.wroc=T){
  x <- substitute(x)
  y <- substitute(y)
  ds <- substitute({
    dat %>%
      group_by(x) %>%
      summarise(goods = sum(y==0),
                bads = sum(y==1))
  }, list(x=x, y=y)) %>%
    eval
  if(output.wroc){
    return(wroc(predictions=ds[[1]],
                labels=ds[c('goods','bads')],
                col.bad = 2, special.values = sp))
  } else{
    return(ds)
  }
}

test_that("wroc's grouped mode can handle various inputs", {
  d <- generate_sample_data(10000, 1000) %>%
    mutate(x_spec2 = round(x_spec),
           x_na2 = round(x_na))

  # test_var <- function
  expect_is(sum_data(d, x, y), 'wroc')
  expect_is(sum_data(d, x_spec2, y, sp=c(-999,-998)), 'wroc')
  expect_is(sum_data(d, x_factor, y), 'wroc')
  expect_is(sum_data(d, x_factor2, y), 'wroc')

  expect_warning(sum_data(d, x_na, y), 'missing|NA|predictions')
  expect_warning(sum_data(d, x_na2, y), 'missing|NA|predictions')
  expect_warning(sum_data(d, x_all_na, y), 'missing|NA|predictions')

  expect_is(sum_data(d, x, y_character), 'wroc')
  expect_is(sum_data(d, x, y_factor), 'wroc')
  expect_warning(sum_data(d, x, y_na), 'NA|missing|labels')
  expect_error(sum_data(d, x, y_all_na), 'labels|missing')

})

test_that("wroc's grouped mode gives the same results as the ungrouped mode", {
  d <- generate_sample_data(10000, 1000) %>%
    mutate(x_spec2 = round(x_spec),
           x_na2 = round(x_na))

  generate_wrocs <- function(dat, x, y, sp=NULL){
    x <- substitute(x)
    y <- substitute(y)
    ds <- substitute({
      dat %>%
        group_by(x) %>%
        summarise(goods = sum(y==0),
                  bads = sum(y==1))
    }, list(x=x, y=y)) %>%
      eval
    suppressWarnings(
      list(
        grouped = wroc(predictions=ds[[1]],
                       labels=ds[c('goods','bads')],
                       col.bad = 2, special.values = sp),
        ungrouped = wroc(predictions=dat[[as.character(x)]],
                         labels=dat[[as.character(y)]],
                         level.bad = 1,
                         special.values = sp)
      )
    )
  }

  expect_equal_wrocs <- function(dat, x, y, sp=NULL){
    x <- substitute(x)
    y <- substitute(y)
    w <- eval(substitute(generate_wrocs(dat, x, y, sp), list(x=x, y=y)))
    expect_equal(w$grouped, w$ungrouped)
  }

  expect_equal_wrocs(d, x, y)
  expect_equal_wrocs(d, x_spec2, y, sp=c(-999,-998))
  expect_equal_wrocs(d, x_na2, y)
  expect_equal_wrocs(d, x_na, y)
  expect_equal_wrocs(d, x_character, y)
  expect_equal_wrocs(d, x_factor2, y)
  expect_equal_wrocs(d, x_all_na, y)
  expect_equal_wrocs(d, x, y_character)
  expect_equal_wrocs(d, x, y_factor)
  expect_equal_wrocs(d, x, y_na)
})

# plot
test_that('plot.wroc works', {
  d <- generate_sample_data(10000, 1000)
  w <- wroc(predictions=d$x_spec, labels=d$y, ngroups = 50,
            special.values = c(-999, -998))
  expect_is(plot(w), 'ggplot')
  expect_is(plot(w, type = 'accum'), 'ggplot')
  expect_is(plot(w, type = 'roc'), 'ggplot')
  expect_is(plot(w, type = 'default'), 'ggplot')
  expect_is(plot(w, type = 'woe'), 'ggplot')

  expect_is(plot(w, type = 'default', include.special = T), 'ggplot')
  expect_is(plot(w, type = 'default', include.special = F), 'ggplot')
  expect_is(plot(w, type = 'woe', include.special = T), 'ggplot')
  expect_is(plot(w, type = 'woe', include.special = F), 'ggplot')
})

# Summary
test_that('summary.wroc works', {
  d <- generate_sample_data(10000, 1000)
  w <- wroc(predictions=d$x_spec, labels=d$y, ngroups = 50,
            special.values = c(-999, -998))
  expect_is(summary(w, performance = T), 'summary.wroc')
  expect_is(summary(w, performance = F), 'summary.wroc')
  expect_null(print(summary(w)))
})


# Analyze

# Optimize



