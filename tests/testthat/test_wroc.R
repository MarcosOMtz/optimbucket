
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

generate_sample_data(10,5)

# Inputs

test_that('wroc can handle all sorts of predictions', {
  d <- generate_sample_data(10000, 1000)

  expect_is(wroc(predictions=d$x, labels=d$y_numeric, ngroups = 50), 'wroc')
  expect_is(wroc(predictions=d$x, labels=d$y_factor, ngroups = 50), 'wroc')
  expect_is(wroc(predictions=d$x, labels=d$y_character, ngroups = 50), 'wroc')
  expect_is(wroc(predictions=d$x, labels=d$y_multilevel, ngroups = 50), 'wroc')
  expect_is(wroc(predictions=d$x, labels=d$y_all_equal, ngroups = 50), 'wroc')

  expect_warning(wroc(predictions=d$x, labels=d$y_na, ngroups = 50), 'NA|labels')
  expect_error(wroc(predictions=d$x, labels=d$y_all_na, ngroups = 50), 'NA|labels')
})

test_that('wroc can handle all sorts of labels', {
  d <- generate_sample_data(10000, 1000)

  expect_is(wroc(predictions=d$x, labels=d$y, ngroups = 50), 'wroc')
  expect_is(wroc(predictions=d$x_factor, labels=d$y, ngroups = 50), 'wroc')
  expect_is(wroc(predictions=d$x_factor2, labels=d$y, ngroups = 50), 'wroc')
  expect_is(wroc(predictions=d$x_character, labels=d$y, ngroups = 50), 'wroc')

  expect_is(wroc(predictions=d$x_spec, labels=d$y, ngroups = 50, special.values = c(-999,-998)), 'wroc')
  expect_warning(wroc(predictions=d$x_na, labels=d$y, ngroups = 50), 'NA|predictions')
  expect_warning(wroc(predictions=d$x_all_na, labels=d$y, ngroups = 50), 'NA|predictions|All|all')
})

# Grouped mode

# Summary

# Analyze

# Optimize



