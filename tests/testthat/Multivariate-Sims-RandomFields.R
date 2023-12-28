
### load other packages
library(MASS)
library(fields)
library(GPvecchia)
library(reshape2)
library(stringr)
library(tidyverse)
library(psych)
library(ggplot2)
#
#
# Multivariate, space-time simulation
# Y(s,t) = X(s) + nu(s,t) + epsilon(s,t)
# Y = Y(i), i = 1,..., k outcomes
#  where:
# X(s) is the spatial only trend,
# nu(s,t) is the spatiotemporally correlated trend, and
# epsilon(s,t) is the spatiotemporal nugget

# To keep the simulation small and tractable to quickly run, we
# simulate 25 spatial locations and 5 temporal locations
n.spatial <- 25
n.spatial.xy <- sqrt(n.spatial)
n.temporal <- 5

# Mean trend, X(s): spatiotemporal data, but the trend is spatial only.
# Each outcome has a different combination of true betas
# p = number of potential covariates
p <- 10

p.nz <- 4
beta1 <- c(rep(0, p - p.nz), rep(1, p.nz))
p.nz <- 2
beta2 <- c(rep(1, p.nz), rep(0, p - p.nz))
p.nz <- 3
beta3 <- c(rep(1, p.nz), rep(0, p - p.nz))


# Combine all of the betas into 1 vector
beta.all <- c(beta1, beta2, beta3)


# correlated predictors
Sigma.X <- exp(-rdist(sample(1:p)) / 3)
X <- mvrnorm(n.spatial, rep(0, p), Sigma.X)
# replicate for space-time for a single-outcomes
X.st <- rbind(X, X, X, X, X)
# Create the multivariate X-matrix, for multiple outcomes
# The rows are by outcome (i), space (j), and then time (k)
X.all <- psych::superMatrix(list(X.st, X.st, X.st))


# S-T design matrix with X's repeated is needed

# Calculate the multivariate  mean trend
# The mean trend is spatial only, but now formatted for space-time (i.e. spatial
# covariates are repeated across time)
mean.trend.st <- X.all %*% beta.all


#### Space-Time correlated error ####
# We simulate using a valid parsimonious Matern spatiotemporal cross-covariance function where
# each marginal parameter (spatial range, temporal range, smoothness, marginal variance)
# are all different. The cross-covariance parameters are the average of the respective
# marginal parameters.
# The simulation is done with unconditional Cholesky decomposition


# rho is the correlation between outcome S-T errors, 1 for with itself
rho <- matrix(nr = 3, c(1, 0.8, 0.3, 0.8, 1, 0.5, 0.3, 0.5, 1))

# nu,  marginal smoothness:
marg.smoothness <- c(0.5, 1, 1.5)

# Non-spatial error/noise/nugget
# epsilon(s,t) will be simulated based on the nugget parameters,
nuggets <- c(0.5, 0.7, 2)

# marginal variances of the Matern
set.seed(10)
# The marginal variances are scaled linearly with the nuggets
x.variance <- runif(3, 1.5, 4)
marg.var <- nuggets + x.variance

# ranges, the true ranges, we will scale coordinates by this, but then use a
# vector of ones in the model/simulation

# outcome1 space and time, outcome 2 space and time, outcome 3 space and time
# Spatial domain is on the unit [0,1]
# Temporal domain is a year [0,365]
ranges <- c(
  0.8, 2,
  0.1, 5,
  0.5, 300
)


# set up the true spatial and temporal dimensions
space.dim <- seq(0, 1, length.out = n.spatial.xy)
time.dim <- seq(1, 365, length.out = n.temporal)

# scale the space and time dimensions according to the true covariance ranges

d11 <- fields::rdist(expand.grid(space.dim / ranges[1], space.dim / ranges[1], time.dim / ranges[2]))
d22 <- fields::rdist(expand.grid(space.dim / ranges[3], space.dim / ranges[3], time.dim / ranges[4]))
d33 <- fields::rdist(expand.grid(space.dim / ranges[5], space.dim / ranges[5], time.dim / ranges[6]))

d12 <- fields::rdist(
  expand.grid(space.dim / ranges[1], space.dim / ranges[1], time.dim / ranges[2]),
  expand.grid(space.dim / ranges[3], space.dim / ranges[3], time.dim / ranges[4])
)
d21 <- t(d12)

d13 <- fields::rdist(
  expand.grid(space.dim / ranges[1], space.dim / ranges[1], time.dim / ranges[2]),
  expand.grid(space.dim / ranges[5], space.dim / ranges[5], time.dim / ranges[6])
)
d31 <- t(d13)

d23 <- fields::rdist(
  expand.grid(space.dim / ranges[3], space.dim / ranges[3], time.dim / ranges[4]),
  expand.grid(space.dim / ranges[5], space.dim / ranges[5], time.dim / ranges[6])
)
d32 <- t(d23)


Dists.All <- rbind(
  cbind(d11, d12, d13),
  cbind(d21, d22, d23),
  cbind(d31, d32, d33)
)

## Create the correlation matrices
Sigma11 <- marg.var[1] * fields::Matern(d11, range = 1, smoothness = marg.smoothness[1])
Sigma22 <- marg.var[2] * fields::Matern(d22, range = 1, smoothness = marg.smoothness[2])
Sigma33 <- marg.var[3] * fields::Matern(d33, range = 1, smoothness = marg.smoothness[3])
Sigma12 <- rho[1, 2] * sqrt(marg.var[1]) * sqrt(marg.var[2]) * fields::Matern(d12, range = 1, smoothness = (marg.smoothness[1] + marg.smoothness[2]) / 2)
Sigma21 <- t(Sigma12)
Sigma13 <- rho[1, 3] * sqrt(marg.var[1]) * sqrt(marg.var[3]) * fields::Matern(d13, range = 1, smoothness = (marg.smoothness[1] + marg.smoothness[3]) / 2)
Sigma31 <- t(Sigma13)
Sigma23 <- rho[2, 3] * sqrt(marg.var[2]) * sqrt(marg.var[3]) * fields::Matern(d23, range = 1, smoothness = (marg.smoothness[2] + marg.smoothness[3]) / 2)
Sigma32 <- t(Sigma23)

# Combine into the super Multivariate covariance matrix
Sigma.All <- rbind(
  cbind(Sigma11, Sigma12, Sigma13),
  cbind(Sigma21, Sigma22, Sigma23),
  cbind(Sigma31, Sigma32, Sigma33)
)

# Cholesky decomposition
L.C <- chol(Sigma.All)
# transpose
L.Ct <- t(L.C)

## Simulate the S-T correlated error
set.seed(10)
st.error <- L.Ct %*% rnorm(n.spatial * n.temporal * 3)

# equivalently, we could do:
# st.error <- rnorm(n.spatial * n.temporal * 3) %*% L.C

# Simulate the nugget (noise) error
set.seed(10)
nug.error <- cbind(
  nuggets[1] * rnorm(n.spatial * n.temporal),
  nuggets[2] * rnorm(n.spatial * n.temporal),
  nuggets[3] * rnorm(n.spatial * n.temporal)
)

# xy and time grid for plotting or accessing by common s-t location
xyt.grid <- cbind(expand.grid(space.dim, space.dim, time.dim))

# outcome indices
idx.outcome1 <- 1:(n.spatial * n.temporal)
idx.outcome2 <- seq(n.spatial * n.temporal + 1, n.spatial * n.temporal * 2, by = 1)
idx.outcome3 <- seq(n.spatial * n.temporal * 2 + 1, n.spatial * n.temporal * 3, by = 1)



###  Combine the mean trend, spatiotemporal errors, and nugget
# Outcome 1, 2, 3 by column

y.sim.final <- cbind(
  mean.trend.st[idx.outcome1] + st.error[idx.outcome1] + nug.error[idx.outcome1],
  mean.trend.st[idx.outcome2] + st.error[idx.outcome2] + nug.error[idx.outcome2],
  mean.trend.st[idx.outcome3] + st.error[idx.outcome3] + nug.error[idx.outcome3]
)

#
###  plot all of the simulations

df.plot <- data.frame(
  "x.coord" = xyt.grid[, 1],
  "y.coord" = xyt.grid[, 2],
  "time" = as.factor(xyt.grid[, 3]),
  "var1" = y.sim.final[, 1],
  "var2" = y.sim.final[, 2],
  "var3" = y.sim.final[, 3],
  "mean.trend1" = mean.trend.st[idx.outcome1],
  "mean.trend2" = mean.trend.st[idx.outcome2],
  "mean.trend3" = mean.trend.st[idx.outcome3],
  "sterror1" = st.error[idx.outcome1],
  "sterror2" = st.error[idx.outcome2],
  "sterror3" = st.error[idx.outcome3],
  "nugerror1" = nug.error[idx.outcome1],
  "nugerror2" = nug.error[idx.outcome2],
  "nugerror3" = nug.error[idx.outcome3]
)


# plot each variable in space and facet by time
ggplot(df.plot, aes(x.coord, y.coord, fill = var1)) +
  geom_raster(interpolate = T) +
  facet_wrap(. ~ time) +
  scale_fill_viridis_c(direction = -1, option = "A")

ggplot(df.plot, aes(x.coord, y.coord, fill = var2)) +
  geom_raster(interpolate = T) +
  facet_wrap(. ~ time) +
  scale_fill_viridis_c(direction = -1, option = "A")

ggplot(df.plot, aes(x.coord, y.coord, fill = var3)) +
  geom_raster(interpolate = T) +
  facet_wrap(. ~ time) +
  scale_fill_viridis_c(direction = -1, option = "A")

# plot the mean trend - each facet should be identical
ggplot(df.plot, aes(x.coord, y.coord, fill = mean.trend1)) +
  geom_raster(interpolate = T) +
  facet_wrap(. ~ time) +
  scale_fill_viridis_c(direction = -1, option = "A")

ggplot(df.plot, aes(x.coord, y.coord, fill = mean.trend2)) +
  geom_raster(interpolate = T) +
  facet_wrap(. ~ time) +
  scale_fill_viridis_c(direction = -1, option = "A")

ggplot(df.plot, aes(x.coord, y.coord, fill = mean.trend3)) +
  geom_raster(interpolate = T) +
  facet_wrap(. ~ time) +
  scale_fill_viridis_c(direction = -1, option = "A")

# plot the st error
ggplot(df.plot, aes(x.coord, y.coord, fill = sterror1)) +
  geom_raster(interpolate = T) +
  facet_wrap(. ~ time) +
  scale_fill_viridis_c(direction = -1, option = "A")

ggplot(df.plot, aes(x.coord, y.coord, fill = sterror2)) +
  geom_raster(interpolate = T) +
  facet_wrap(. ~ time) +
  scale_fill_viridis_c(direction = -1, option = "A")

# Note how the S-T error 3 is the smoothest and also the most similar across time.
# It has the largest smoothness and temporal range parameters
ggplot(df.plot, aes(x.coord, y.coord, fill = sterror3)) +
  geom_raster(interpolate = T) +
  facet_wrap(. ~ time) +
  scale_fill_viridis_c(direction = -1, option = "A")

# plot the nugget
ggplot(df.plot, aes(x.coord, y.coord, fill = nugerror1)) +
  geom_raster(interpolate = T) +
  facet_wrap(. ~ time) +
  scale_fill_viridis_c(direction = -1, option = "A")

ggplot(df.plot, aes(x.coord, y.coord, fill = nugerror2)) +
  geom_raster(interpolate = T) +
  facet_wrap(. ~ time) +
  scale_fill_viridis_c(direction = -1, option = "A")

ggplot(df.plot, aes(x.coord, y.coord, fill = nugerror3)) +
  geom_raster(interpolate = T) +
  facet_wrap(. ~ time) +
  scale_fill_viridis_c(direction = -1, option = "A")


# Some simulation field stats

# Signal to noise
S2N.var1 <- var(df.plot$mean.trend1) / var(df.plot$sterror1 + df.plot$nugerror1)
S2N.var2 <- var(df.plot$mean.trend2) / var(df.plot$sterror2 + df.plot$nugerror2)
S2N.var3 <- var(df.plot$mean.trend3) / var(df.plot$sterror3 + df.plot$nugerror3)

X_train <- as.matrix(X.st[1:90, ])
X_test <- as.matrix(X.st[91:125, ])
Y_train <- as.matrix(y.sim.final[1:90, ])
Y_test <- as.matrix(y.sim.final[91:125, ])
# locs_train = xyt.grid[1:90,]
# locs_test = xyt.grid[91:125,]
noise_sd <- sd(xyt.grid[, 1]) * 0.01
noise <- cbind(rnorm(nrow(xyt.grid), 0, noise_sd), rnorm(nrow(xyt.grid), 0, noise_sd), rnorm(nrow(xyt.grid), 0, sd(xyt.grid[, 3]) * 0.01))
xyt_jittered <- xyt.grid + noise
locs_train <- as.matrix(xyt_jittered[1:90, 1:2]) # spatial only
locs_test <- as.matrix(xyt_jittered[91:125, 1:2]) # spatial only
save(X_train, X_test, Y_train, Y_test, locs_train, locs_test, file = "tests/testthat/multivariate_sim_spatial3.Rdata")
