###  a simple simulated example for variable section + general vecchia


### load other packages
library(MASS)
library(fields)




##### true parameters
n <- 2000
p <- 50
p.nz <- 5 # Number of true non-zero coefficients
beta0 <- c(rep(1, p.nz), rep(0, p - p.nz))
theta0 <- c(9, .3, .5, 1) # variance,range,smoothness,nugget


##### simulate data
locs <- cbind(runif(n), runif(n)) # random sample on unit square
Sigma.X <- exp(-rdist(sample(1:p)) / 5)
X <- mvrnorm(n, rep(0, p), Sigma.X) # correlated predictors

# Error terms for simulation that includes spatial and nugget
Sigma0 <- theta0[1] * Matern(rdist(locs), range = theta0[2], smoothness = theta0[3]) + theta0[4] * diag(n)
epsilon <- mvrnorm(1, mu = rep(0, n), Sigma0) # Simulate the GP
y <- X %*% beta0 + epsilon # Calculate the observed data

# If we want the true underlying field, then we could simply do X %*% beta0

X_train <- X[1:1000, ]
X_test <- X[1001:2000, ]
Y_train <- y[1:1000, ]
Y_test <- y[1001:2000, ]
locs_train <- locs[1:1000, ]
locs_test <- locs[1001:2000, ]
save(X_train, X_test, Y_train, Y_test, locs_train, locs_test, file = "tests/testthat/sim_spatial.Rdata")
