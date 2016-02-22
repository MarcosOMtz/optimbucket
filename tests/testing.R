#### Datos 1
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

#### Datos 2
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

#### Datos 3
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

#### Ejemplo

wr <- wroc(d$x, d$y, ngroups=20)
auc.wroc(wr)
plot(wr)
plot(wr, type='trend')
plot(wr, type='woe')

wr2 <- optimize.wroc(wr, 'auto')
auc.wroc(wr2)
plot(wr2)
plot(wr2, type='trend')
plot(wr2, type='woe')

# Trameado manual usando subset
wr3 <- subset(wr, c(1:5,7:17,19:26,28:33,35:36,38:39,43))
plot(wr3, type='woe')
