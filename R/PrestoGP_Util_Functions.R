################################################################################
# SCAD Penalty value
################################################################################
SCAD_Penalty_Loglike <- function(beta.in, lambda) {
  penalty <- matrix(NA, nrow = nrow(beta.in), ncol = ncol(beta.in))
  for (j in 1:nrow(beta.in)) {
    beta.j <- beta.in[j, ]
    idx1 <- abs(beta.j) <= lambda
    idx2 <- abs(beta.j) > lambda & abs(beta.j) <= 3.7 * lambda
    idx3 <- abs(beta.j) > 3.7 * lambda
    penalty[j, idx1] <- lambda[idx1] * beta.j[idx1]
    penalty[j, idx2] <- -(abs(beta.j[idx2])^2 - 7.4 * lambda[idx2] * abs(beta.j[idx2]) + lambda[idx2]^2) / (5.4)
    penalty[j, idx3] <- (3.7 * lambda[idx3]^2 + lambda[idx3]^2) / 2
  }


  loglik.penalty <- lambda * colSums(penalty)

  return(loglik.penalty)
}

#' glmnet_penalty
#'
#' Compute a combination of L1 and L2 penalties
#'
#' @param beta the beta coefficients
#' @param lambda the penalty factor (how big the total penalty is)
#' @param alpha the L1 vs L2 coefficient (0 means all L2, 1 means all L1)
#'
#' @return a penalty score
glmnet_penalty <- function(beta, lambda, alpha) {
  return(lambda * ((1 - alpha) * sqrt(sum(beta^2)) + alpha * sum(abs(beta))))
}

#' Ordinary Spatiotemporal Kriging Prediction with a Local S-T neigborhood.
#'
#' @param new_coords Locations to predict
#' @param obs_coords Training locations
#' @param Y_obs Training dependent variable vector
#' @param cov.pars Estimated covariance parameters
#' @param NN number of nearest neighbors to use
#'
#' @return A dataframe with predicted means and variances
#' @export
#'
#' @examples
#' @noRd
Kr_pred <- function(new_coords, obs_coords, Y_obs, cov.pars, NN) {
  Kr.prediction <- matrix(NA, nrow = nrow(new_coords), ncol = 1)
  Kr.Var <- matrix(NA, nrow = nrow(new_coords), ncol = 1)

  new_coords_scaled <- cbind(new_coords[, 1] / cov.pars[2], new_coords[, 2] / cov.pars[2], new_coords[, 3] / cov.pars[3])
  obs_coords_scaled <- cbind(obs_coords[, 1] / cov.pars[2], obs_coords[, 2] / cov.pars[2], obs_coords[, 3] / cov.pars[3])
  df_new <- new_coords_scaled
  df_obs <- obs_coords_scaled

  for (i in 1:nrow(new_coords)) {
    # print(i)
    locs.test.i <- rep.row(df_new[i, ], nrow(df_obs))
    dist.op <- fields::rdist.vec(locs.test.i, df_obs)
    Sigma.op <- cov.pars[1] * fields::Exponential(dist.op, range = 1)
    s.op <- sort(dist.op, index.return = TRUE)
    s.idx <- s.op$ix[1:NN]
    # Get the closest "NN" observations
    dist.oo <- fields::rdist(df_obs[s.idx, ], df_obs[s.idx, ])
    Sigma.oo <- cov.pars[4] * diag(NN) +
      cov.pars[1] * fields::Exponential(dist.oo, range = 1)
    oo <- Y_obs[s.idx]
    ## kriging predictor (posterior mean)
    mean_trend <- mean(oo)
    Kr.prediction[i] <- mean_trend + t(Sigma.op[s.idx]) %*% solve(Sigma.oo, oo - mean_trend)

    Kr.Var[i] <- cov.pars[1] - t(Sigma.op[s.idx]) %*% solve(Sigma.oo, Sigma.op[s.idx])
  }

  Kr.Data <- data.frame(Kr.prediction, Kr.Var)

  return(Kr.Data)
}

#' ST_Krig_Param_Avg
#'
#' Spatiotemporal Kriging Maximum Likelihood Estiamtion based on an average of multiple random subsets.
#'
#' @param Y A vector containing training values for the dependent variable.
#' @param locs A matrix containing the training spatial coordinates and times.
#' @param p The number of data points to sample with each iteration.
#' @param k The number of optimization iterations.
#'
#' @return a vector containing the covariance parameters
#' @export
#'
#' @examples
#' @noRd
ST_Krig_Param_Avg <- function(Y, locs, p, k = 10) {
  n <- length(Y)
  mdl.geo.fit.avg <- list()
  mdl.geo.fit.avg <- matrix(NA, nrow = k, ncol = 4)

  for (i in 1:k) {
    set.seed(i)
    OK.fold <- sample(1:n, p, replace = FALSE)
    Y.train <- Y[OK.fold]
    locs.train <- locs[OK.fold, ]
    D.sample <- rdist(locs.train[, 1:2])
    t.sample <- rdist(locs.train[, 3])
    theta.hat <- c(.9 * var(Y.train), mean(D.sample) / 4, mean(t.sample) / 4, 0.1 * var(Y.train)) # var,s-range,t-range,nugget
    full.result <- optim(
      par = log(theta.hat), fn = negloglik_full_ST,
      locs = locs.train, y = Y.train, N = p,
      control = list(trace = FALSE, maxit = 400)
    )
    mdl.geo.fit.avg[i, ] <- full.result$par
  }

  covparam <- exp(apply(mdl.geo.fit.avg, 2, mean))

  return(covparam)
}


################################################################################
## replicate rows, helps with vector distance
################################################################################
rep.row <- function(x, n) {
  matrix(rep(x, each = n), nrow = n)
}


################################################################################
## Transform Spatiotemporal data to i.i.d.
################################################################################
transform_iid <- function(data, vecchia.approx, covparms, nuggets) {
  # compute required matrices
  U.obj <- createU(vecchia.approx, covparms, nuggets)
  V.ord <- U2V(U.obj)
  U.z <- U.obj$U[!U.obj$latent, ]
  U.y <- U.obj$U[U.obj$latent, ]

  # compute transformed data in parts
  part1.ord <- Matrix::crossprod(U.z, data[U.obj$ord.z, ]) # C.hat^-1
  temp1 <- U.y %*% part1.ord #
  revord <- nrow(temp1):1
  temp2 <- spam::solve(V.ord, temp1[revord, ])
  part2.rev <- spam::solve(Matrix::t(V.ord), temp2)
  part2.ord <- spam::crossprod(U.y, part2.rev[revord, ])
  transform.ord <- part1.ord - part2.ord

  # return to original ordering
  orig.order <- order(U.obj$ord)
  transformed.data <- transform.ord[orig.order, ]
  return(transform.ord)
}

################################################################################
## Transform Multivariate Spatiotemporal data to i.i.d.
################################################################################
transform_miid <- function(data, vecchia.approx, params) {
  # compute required matrices
  U.obj <- createUMultivariate(vecchia.approx, params)
  if (sum(is.na(U.obj$U)) > 0) {
    browser()
  }
  V.ord <- U2V(U.obj)
  U.z <- U.obj$U[!U.obj$latent, ]
  U.y <- U.obj$U[U.obj$latent, ]

  # compute transformed data in parts
  part1.ord <- Matrix::crossprod(U.z, data[U.obj$ord.z, ]) # C.hat^-1
  temp1 <- U.y %*% part1.ord #
  revord <- nrow(temp1):1
  temp2 <- spam::solve(V.ord, temp1[revord, ])
  part2.rev <- spam::solve(Matrix::t(V.ord), temp2)
  part2.ord <- spam::crossprod(U.y, part2.rev[revord, ])
  transform.ord <- part1.ord - part2.ord

  # return to original ordering
  orig.order <- order(U.obj$ord)
  transformed.data <- transform.ord[orig.order, ]
  return(transform.ord)
}

######  GPvecchia local function
###### compute V for posterior inference - needed for transform.iid   #######
U2V <- function(U.obj) {
  U.y <- U.obj$U[U.obj$latent, ]

  if (U.obj$cond.yz == "zy") {
    V.ord <- revMat(U.y[, U.obj$latent, drop = FALSE])
  } else if (U.obj$ord.pred != "obspred") {
    W <- Matrix::tcrossprod(U.y)
    W.rev <- revMat(W)
    V.ord <- Matrix::t(Matrix::chol(W.rev))
  } else { # for obspred ordering

    last.obs <- max(which(!U.obj$latent))
    latents.before <- sum(U.obj$latent[1:last.obs])
    latents.after <- sum(U.obj$latent[-(1:last.obs)])

    # pred columns are unchanged
    V.pr <- revMat(U.y[, (last.obs + 1):ncol(U.y), drop = FALSE])

    # have to compute cholesky for obs block
    U.oo <- U.y[1:latents.before, 1:last.obs]
    A <- Matrix::tcrossprod(U.oo)
    A.rev <- revMat(A)
    V.oor <- Matrix::t(Matrix::chol(A.rev))

    # combine the blocks into one matrix
    zeromat.sparse <- Matrix::sparseMatrix(c(), c(), dims = c(latents.after, latents.before))
    V.or <- rbind(zeromat.sparse, V.oor)

    V.ord <- methods::as(cbind(V.pr, V.or), "triangularMatrix")
  }

  return(V.ord)
}

###########################################################
## Reverse order of matrix rows,cols
revMat <- function(mat) mat[nrow(mat):1, ncol(mat):1, drop = FALSE]

#' @export
vecchia_Mprediction <- function(z, vecchia.approx, covparms, var.exact = NULL,
                                return.values = "mean") {
  GPvecchia:::removeNAs()
  U.obj <- createUMultivariate(vecchia.approx, covparms)
  V.ord <- U2V(U.obj)
  #    if (length(U.obj$zero.nugg) > 0)
  #        warning("Rows/cols of V have been removed for data with zero noise")
  vecchia.mean <- GPvecchia:::vecchia_mean(z, U.obj, V.ord)
  return.list <- list(
    mu.pred = vecchia.mean$mu.pred, mu.obs = vecchia.mean$mu.obs,
    var.pred = NULL, var.obs = NULL, V.ord = NULL, U.obj = NULL
  )
  if (return.values == "meanmat" | return.values == "all") {
    return.list$V.ord <- V.ord
    return.list$U.obj <- U.obj
  }
  if (return.values == "meanvar" | return.values == "all") {
    if (is.null(var.exact)) {
      var.exact <- (sum(!vecchia.approx$obs) < 4 * 10000)
    }
    vars.vecchia <- GPvecchia:::vecchia_var(U.obj, V.ord, exact = var.exact)
    return.list$var.pred <- vars.vecchia$vars.pred
    return.list$var.obs <- vars.vecchia$vars.obs
  }
  return(return.list)
}
