set.seed(1212)

ny <- 3 # number of response variables
n.spatial.xy <- 10 # number of spatial coordinates per dimension

library(MASS)
library(fields)
library(psych)

n.rho <- choose(ny, 2)
rho.vec <- runif(n.rho, 0.2, 0.8)
rho <- matrix(0, nrow = ny, ncol = ny)
rho[upper.tri(rho)] <- rho.vec
rho <- rho + t(rho) + diag(1, ny)

marg.smoothness <- 0.5 + rnorm(ny, 0.1, 0.05)
nuggets <- runif(ny, 0.5, 2)
x.variance <- runif(ny, 1.5, 4)
marg.var <- nuggets + x.variance
ranges <- runif(ny, 0.5, 1.2)

params.all <- c(x.variance, ranges, marg.smoothness, nuggets, rho.vec)

locs.list <- list()
for (i in 1:ny) {
  loc1 <- seq(0, 1, length.out = n.spatial.xy) + rnorm(n.spatial.xy, 0, 0.001)
  loc2 <- seq(0, 1, length.out = n.spatial.xy) + rnorm(n.spatial.xy, 0, 0.001)
  locs.list[[i]] <- as.matrix(expand.grid(loc1, loc2))
}

Sigma.All <- matrix(nrow = ny * n.spatial.xy^2, ncol = ny * n.spatial.xy^2)
for (i in 1:ny) {
  for (j in i:ny) {
    if (i == j) {
      ndx1 <- ((i - 1) * n.spatial.xy^2 + 1):(n.spatial.xy^2 * i)
      dij <- fields::rdist(locs.list[[i]])
      Sigma.All[ndx1, ndx1] <- marg.var[i] * fields::Matern(dij, range = ranges[i], smoothness = marg.smoothness[i])
    } else {
      ndx1 <- ((i - 1) * n.spatial.xy^2 + 1):(n.spatial.xy^2 * i)
      ndx2 <- ((j - 1) * n.spatial.xy^2 + 1):(n.spatial.xy^2 * j)
      dij <- fields::rdist(locs.list[[i]], locs.list[[j]])
      vii <- marg.smoothness[i]
      vjj <- marg.smoothness[j]
      vij <- (vii + vjj) / 2
      aii <- 1 / ranges[i]
      ajj <- 1 / ranges[j]
      aij <- sqrt((aii^2 + ajj^2) / 2)
      Sigma.All[ndx1, ndx2] <- rho[i, j] * sqrt(marg.var[i]) * sqrt(marg.var[j]) * aii^vii * ajj^vjj * gamma(vij) /
        (aij^(2 * vij) * sqrt(gamma(vii) * gamma(vjj))) * fields::Matern(dij, smoothness = vij, alpha = aij)
      Sigma.All[ndx2, ndx1] <- t(Sigma.All[ndx1, ndx2])
    }
  }
}

L.C <- chol(Sigma.All)

st.error <- rnorm(n.spatial.xy^2 * ny) %*% L.C

nug.error <- NULL
for (i in 1:ny) {
  nug.error <- c(nug.error, nuggets[i] * rnorm(n.spatial.xy^2))
}

y.list <- list()
for (i in 1:ny) {
  ndx1 <- ((i - 1) * n.spatial.xy^2 + 1):(n.spatial.xy^2 * i)
  y.list[[i]] <- st.error[ndx1] + nug.error[ndx1]
}

rm(
  ny, n.spatial.xy, i, j, n.rho,
  loc1, loc2, ndx1, ndx2, dij, vii, vjj, vij, aii, ajj, aij, L.C,
  st.error, nug.error, rho, rho.vec, ranges, Sigma.All, nuggets,
  marg.smoothness, marg.var
)
