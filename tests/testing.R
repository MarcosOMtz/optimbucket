require(dplyr)
require(tidyr)
require(ggplot2)

#### Data 1
ng <- 10000
nb <- 1000

goods <- data.frame(
  x = round(rnorm(ng, mean = -2),2),
  y = 0
)
bads <- data.frame(
  x = rnorm(nb),
  y = 1
)

d <- rbind(goods, bads) %>%
  mutate(y = factor(y))

ggplot(d, aes(x, fill=y)) +
  geom_density(alpha=0.5)

#### Data 2
ng <- 1000
nb <- 100

goods <- data.frame(
  x = rnorm(ng, mean = ifelse(runif(ng) > 0.8, 4, -2)),
  y = 0
)
bads <- data.frame(
  x = rnorm(nb),
  y = 1
)

d <- rbind(goods, bads) %>%
  mutate(y = factor(y))

ggplot(d, aes(x, fill=y)) +
  geom_density(alpha=0.5)

#### Data 3
ng <- 10000
nb <- 1000

goods <- data.frame(
  x = rnorm(ng, mean = ifelse(runif(ng) > 0.2, 10, -10)),
  y = 0
)
bads <- data.frame(
  x = rnorm(nb),
  y = 1
)

d <- rbind(goods, bads) %>%
  mutate(y = factor(y))

ggplot(d, aes(x, fill=y)) +
  geom_density(alpha=0.5)

#### Data 4
ng <- 10000
nb <- 1000

goods <- data.frame(
  x = rnorm(ng, mean = ifelse(runif(ng) > 0.2, -10, 10)),
  y = 0
)
bads <- data.frame(
  x = rnorm(nb),
  y = 1
)

d <- rbind(goods, bads) %>%
  mutate(y = factor(y))

ggplot(d, aes(x, fill=y)) +
  geom_density(alpha=0.5)

#### Data 5
ng <- 10000
nb <- 1000

goods <- data.frame(
  x = rnorm(ng, mean = ifelse(runif(ng) > 0.2, -2, 1)),
  y = 0
)
bads <- data.frame(
  x = rnorm(nb),
  y = 1
)

d <- rbind(goods, bads) %>%
  mutate(y = factor(y))

ggplot(d, aes(x, fill=y)) +
  geom_density(alpha=0.5)

#### Data 6 MULTI
ng <- 10000
nb <- 1000

goods <- data.frame(
  x = rnorm(ng, mean = ifelse(runif(ng) > 0.2, -2, 1)),
  y = 0
)
bads <- data.frame(
  x = rnorm(nb),
  y = 1
)

d <- rbind(goods, bads) %>%
  mutate(y = factor(y),
         z = 2*x + 1 + rnorm(ng+nb, 0, 5),
         w = round(x - z  - 1 + rnorm(ng+nb, 0, 1)),
         w = ifelse(runif(ng+nb) < 0.01,
                    -999,
                    ifelse(runif(ng+nb) < 0.01,
                           -998,
                           w)))

ggplot(d, aes(x, fill=y)) +
  geom_density(alpha=0.5)

ggplot(d, aes(z, fill=y)) +
  geom_density(alpha=0.5)

ggplot(d, aes(w, fill=y)) +
  geom_density(alpha=0.5)

#### Example

# Create and examine a  wroc object
wr <- wroc(d$x, d$y, ngroups=20)
summary(wr) # Just the global details
summary(wr, performance = F) # Faster but less informative
print(summary(wr), extended = TRUE)
performance(wr)
plot(wr)
plot(wr, type='trend')
plot(wr, type='woe')

# Optimal Gini bucketing
wr2 <- optimize(wr, 'auto')
performance(wr2)
plot(wr2)
plot(wr2, type='trend')
plot(wr2, type='woe')

# Manual bucketing
wr3 <- subset(wr, c(1:5,7:17))
plot(wr3, type='woe')

# Analysis of the optimal path
plot(wr4 <- analyze(wr, 7), 'woe')

# wroc's for multiple variables at once
wrs <- wroc(y ~ x + z, d, ngroups = 50, level.bad = 1,
            special.values = list(w = c(-998, -999)))

# wroc for all variables
wrs <- wroc(y ~ ., d, ngroups = 50, level.bad = 1,
            special.values = list(w = c(-998, -999)))
plot(wrs, save.pdf = F)

# Manipulate wroc.list objects
subset(wrs, keep=NULL, drop='z')
wrs[c('x','z')]
wrs['z']
wrs[2]

wrs2 <- wrs
names(wrs2) <- c('a','b','c')
wrs3 <- c(wrs, wrs2)

# Optimize a bunch of variables
owrs <- optimize(wrs, 'auto')
owrs

# Use a wroc object to paste WoE (or any other variable in a wroc$info table, such as bucket number or probability of default)
predict(object = wrs$w, newdata = head(d),
        variable = 'w', type = 'woe', keep.data = T)

# Use a wroc.list object to paste several WoEs in one step
predict(object = owrs, newdata = head(d), type = 'woe', keep.data = T)











