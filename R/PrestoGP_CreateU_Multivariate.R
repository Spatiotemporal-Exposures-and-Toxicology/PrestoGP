##############################################################################
### Create the parameter sequence indices for the likelihood function #########

create.param.sequence <- function(P, nr = 1) {
  # Input: P - number of multivariate outcome dimensions
  # Assumes our version of the flexible multivariate matern model
  nk <- choose(P, 2)
  if (nk == 0) {
    nk <- 1 # univariate case
  }

  #    param.sequence.begin <- c(seq(1,P*4,by=P),seq(P*4+1,(P*4)+(nk),by=nk))
  param.sequence.begin <- c(1, P + 1, seq(P * (nr + 1) + 1, length = 3, by = P))
  param.sequence.end <- c(P, nr * P, P, P, nk) %>% cumsum()
  param.sequence <- cbind(param.sequence.begin, param.sequence.end)


  return(param.sequence)
}

#' max_min_ordering
#'
#' Determine an ordering of locations that uses the max-min method to get a sparse covering of the location space.
#'
#' @param locs A matrix with one row per location and any number of columns (x, y, time, etc).
#' @param dist_func Any distance function with a signature of dist(query_location, locations_matrix)
#'
#' @return A vector containing the index of locations from sparse to dense
max_min_ordering <- function(locs, dist_func) {
  center <- matrix(colMeans(locs), ncol = ncol(locs))
  # find the point closest to the mean of all points
  dists <- dist_func(center, locs)
  first <- which.min(dists)
  unsolved <- 1:nrow(locs)
  unsolved <- unsolved[-first]
  order <- c(first)

  while (length(order) < nrow(locs)) {
    max_min <- 0
    max_min_i <- unsolved[1]
    in_order <- locs[order[1:length(order)], ]
    dim(in_order) <- c(length(order), ncol(locs))
    for (i in unsolved) {
      loc_i <- locs[i, ]
      dim(loc_i) <- c(1, ncol(locs))
      dists <- dist_func(loc_i, in_order)
      candiate_dist <- min(dists) # min distance from loc(i) to already ordered points
      if (candiate_dist > max_min) {
        max_min <- candiate_dist # loc(i) has a larger minimum distance from already ordered points
        max_min_i <- i
      }
    }
    order <- c(order, max_min_i)
    unsolved <- unsolved[-which(unsolved == max_min_i)] # mark the max-min loc as solved
  }
  order
}

#' knn_indices
#'
#' Find the index of K nearest neighbors within a set of locations for a given query and distance function
#'
#' @param ordered_locs A matrix with one row per location, where locations are ordered using max_min ordering.
#' @param query Find the nearest neighbors to this location
#' @param n_neighbors The number of neighbors to find (K)
#' @param dist_func Any distance function with a signature of dist(query_location, locations_matrix)
#'
#' @return A vector containing the indices of the neighbors
knn_indices <- function(ordered_locs, query, n_neighbors, dist_func, dist_func_code) {
  if (dist_func_code == "custom") {
    dists <- dist_func(query, ordered_locs)
    dists_order <- order(dists)
    nearest_neighbors <- dists_order[1:n_neighbors]
    return(list(
      "indices" = nearest_neighbors,
      "distances" = dists[nearest_neighbors]
    ))
  } else {
    cur.nn <- RANN::nn2(ordered_locs, query, n_neighbors)
    return(list("indices" = cur.nn$nn.idx, "distances" = cur.nn$nn.dists))
  }
}

#' sparseNN
#'
#' Find the index of K nearest neighbors within a set of locations for a given query and distance function
#'
#' @param ordered_locs A matrix with one row per location, where locations are ordered using max_min ordering.
#' @param n_neighbors The number of neighbors to find (K) for each location
#' @param dist_func Any distance function with a signature of dist(query_location, locations_matrix)
#'
#' @return A list containing two matrices, each with one row per location: an indices matrix with the indices of nearest neighbors for each location, and a distance matrix with the associated distances
sparseNN <- function(ordered_locs, n_neighbors, dist_func, dist_func_code, ordered_locs_pred = NULL) {
  ee <- min(apply(ordered_locs, 2, stats::sd))
  n <- nrow(ordered_locs)
  ordered_locs <- ordered_locs + matrix(
    ee * 1e-04 *
      stats::rnorm(n * ncol(ordered_locs)),
    n, ncol(ordered_locs)
  )
  indices_matrix <- matrix(data = NA, nrow = nrow(ordered_locs), ncol = n_neighbors)
  distances_matrix <- matrix(data = NA, nrow = nrow(ordered_locs), ncol = n_neighbors)
  for (row in 1:n_neighbors) {
    # for the locations from 1 to n_neighbors, use the entire locs list to find the neighbors
    nn <- knn_indices(ordered_locs[1:(n_neighbors + 1), , drop = FALSE][-row, , drop = FALSE], ordered_locs[row, , drop = FALSE], n_neighbors, dist_func, dist_func_code)
    indices_matrix[row, 1:n_neighbors] <- nn$indices[1:n_neighbors]
    distances_matrix[row, 1:n_neighbors] <- nn$distances[1:n_neighbors]
  }
  for (row in (n_neighbors + 1):nrow(ordered_locs)) {
    # get the m nearest neighbors from the locs before this one in the max-min order
    nn <- knn_indices(ordered_locs[1:(row - 1), , drop = FALSE], ordered_locs[row, , drop = FALSE], n_neighbors, dist_func, dist_func_code)
    indices_matrix[row, 1:n_neighbors] <- nn$indices[1:n_neighbors]
    distances_matrix[row, 1:n_neighbors] <- nn$distances[1:n_neighbors]
  }
  if (!is.null(ordered_locs_pred)) {
    indices_matrix_pred <- matrix(
      data = NA, nrow = nrow(ordered_locs_pred),
      ncol = n_neighbors
    )
    distances_matrix_pred <- matrix(
      data = NA, nrow = nrow(ordered_locs_pred),
      ncol = n_neighbors
    )
    for (row in 1:nrow(ordered_locs_pred)) {
      nn <- knn_indices(
        ordered_locs,
        ordered_locs_pred[row, , drop = FALSE], n_neighbors,
        dist_func, dist_func_code
      )
      indices_matrix_pred[row, 1:n_neighbors] <- nn$indices[1:n_neighbors]
      distances_matrix_pred[row, 1:n_neighbors] <- nn$distances[1:n_neighbors]
    }
    indices_matrix <- rbind(indices_matrix, indices_matrix_pred)
    distances_matrix <- rbind(distances_matrix, distances_matrix_pred)
  }
  list("indices" = indices_matrix, "distances" = distances_matrix)
}

# This partitions the vector q(i) of nearest neighbors into the subsets
# q_y(i) and q_z(i) as described in the second full paragraph on page
# 17 of Kyle's technical report.
calc.q <- function(nn.obj, firstind.pred) {
  m <- ncol(nn.obj)
  n <- nrow(nn.obj)

  q.y <- list(length = n)
  q.z <- list(length = n)
  q.y[[1]] <- NULL
  q.z[[1]] <- NULL

  for (i in 2:(m + 1)) {
    q.y[[i]] <- 1:(i - 1)
    q.z[[i]] <- NULL
  }

  for (i in (m + 1):n) {
    cur.q <- nn.obj[i, ]
    best.k <- cur.q[1]
    best.qy <- intersect(q.y[[best.k]], cur.q)
    for (j in 2:m) {
      cur.k <- cur.q[j]
      cur.qy <- intersect(q.y[[cur.k]], cur.q)
      if (length(cur.qy) > length(best.qy) & cur.k < firstind.pred) {
        best.k <- cur.k
        best.qy <- cur.qy
      }
    }
    q.y[[i]] <- union(best.k, best.qy)
    q.z[[i]] <- setdiff(cur.q, q.y[[i]])
    latent.pred <- q.z[[i]][q.z[[i]] >= firstind.pred]
    if (length(latent.pred) > 0) {
      q.y[[i]] <- union(q.y[[i]], latent.pred)
      q.z[[i]] <- setdiff(q.z[[i]], latent.pred)
    }
  }
  return(list(q.y = q.y, q.z = q.z))
}

#' @export
vecchia_Mspecify <- function(locs.list, m, locs.list.pred = NULL,
                             dist.func = NULL,
                             ordering.pred = c("obspred", "general"),
                             pred.cond = c("independent", "general")) {
  ordering.pred <- match.arg(ordering.pred)
  pred.cond <- match.arg(pred.cond)

  dist.func.code <- "custom"
  if (is.null(dist.func)) {
    dist.func <- fields::rdist
    dist.func.code <- "rdist"
  }

  P <- length(locs.list)

  locs <- NULL
  ndx <- NULL
  for (i in 1:P) {
    if (!is.matrix(locs.list[[i]])) {
      stop("Each element of locs.list must be a matrix")
    }
    locs <- rbind(locs, locs.list[[i]])
    ndx <- c(ndx, rep(i, nrow(locs.list[[i]])))
  }
  n <- nrow(locs)

  if (!is.null(locs.list.pred)) {
    if (length(locs.list.pred) != P) {
      stop("locs.list and locs.list.pred must have the same length")
    }
    locs.pred <- NULL
    ndx.pred <- NULL
    for (i in 1:P) {
      if (!is.null(locs.list.pred[[i]])) {
        if (!is.matrix(locs.list.pred[[i]])) {
          stop("Each element of locs.list.pred must be a matrix or NULL")
        }
        locs.pred <- rbind(locs.pred, locs.list.pred[[i]])
        ndx.pred <- c(ndx.pred, rep(i, nrow(locs.list.pred[[i]])))
      }
    }
    locs.all <- rbind(locs, locs.pred)
    ndx.all <- c(ndx, ndx.pred)
    observed.obspred <- c(rep(TRUE, n), rep(FALSE, nrow(locs.pred)))
  } else {
    locs.all <- locs
    ndx.all <- ndx
    observed.obspred <- rep(TRUE, n)
    ordering.pred <- "general"
  }

  if (dist.func.code == "custom") {
    if (!is.null(locs.list.pred)) {
      stop("Only Euclidean distance currently supported for prediction")
    }
    loc.order <- max_min_ordering(locs.all, dist.func)
    loc.order <- c(unique(loc.order), setdiff(1:n, loc.order))
  } else {
    if (is.null(locs.list.pred) | ordering.pred == "general") {
      loc.order <- GPvecchia::order_maxmin_exact(locs.all)
      # I am not sure why the next two lines are here. I added them because
      # similar code exists in the GPvecchia package. But I don't know why
      # they did this.
      cutoff <- min(n, 9)
      loc.order <- c(
        loc.order[1], loc.order[-seq(1, cutoff)],
        loc.order[2:cutoff]
      )
      ord <- loc.order
      ord.z <- loc.order[loc.order <= n]
    } else {
      loc.order <- GPvecchia::order_maxmin_exact_obs_pred(
        locs,
        locs.pred
      )
      ord.z <- loc.order$ord
      ord <- c(ord.z, loc.order$ord_pred + nrow(locs))
    }
    olocs <- locs.all[ord, , drop = FALSE]
    ondx <- ndx.all[ord]
    obs <- observed.obspred[ord]
  }

  # Note that the corresponding function in the GPvecchia package
  # is non-deterministic, so there may be some slight differences
  # between the output of this function and the output of createU
  # in the GPvecchia package.
  if (is.null(locs.list.pred) | pred.cond == "general") {
    nn.mat <- sparseNN(olocs, m, dist.func, dist.func.code)
  } else {
    nn.mat <- sparseNN(
      olocs[1:n, , drop = FALSE], m, dist.func, dist.func.code,
      olocs[-(1:n), , drop = FALSE]
    )
  }
  last.obs <- max((1:length(obs))[obs])
  q.list <- calc.q(nn.mat$indices, last.obs + 1)

  return(list(
    locsord = olocs, obs = obs, ord = ord, ord.z = ord.z,
    ord.pred = ordering.pred, cond.yz = "SGV", conditioning = "NN",
    P = P, ondx = ondx, dist.func = dist.func,
    dist.func.code = dist.func.code, q.list = q.list,
    n.neighbors = m
  ))
}

# This function computes the U matrix for a Vecchia approximation similar
# to the createU function from the GPvecchia package. The critical
# difference is that it can handle multivariate input data.
#
# locs.list is a list of coordinates (one for each outcome variable)
# params is a vector of parameters that is split into its relevant
# parts using the create.param.sequence function. Note that currently
# it allows only one scale parameter per outcome. This needs to be
# fixed.
#' @useDynLib PrestoGP
#' @export
createUMultivariate <- function(vec.approx, params, cov_func = NULL) {
  if (is.null(cov_func)) {
    cov_func <- fields::Matern
  }

  dist_func <- vec.approx$dist.func
  P <- vec.approx$P
  q.list <- vec.approx$q.list
  olocs <- vec.approx$locsord
  n <- nrow(olocs)
  ondx <- vec.approx$ondx

  param.seq <- create.param.sequence(P)
  param.sequence.begin <- param.seq[, 1]
  param.sequence.end <- param.seq[, 2]

  sig2 <- params[param.sequence.begin[1]:param.sequence.end[1]]
  rangep <- params[param.sequence.begin[2]:param.sequence.end[2]]
  smoothness <- params[param.sequence.begin[3]:param.sequence.end[3]]
  nugget <- params[param.sequence.begin[4]:param.sequence.end[4]]

  rho <- params[param.sequence.begin[5]:param.sequence.end[5]]
  rho.mat <- matrix(0, nrow = P, ncol = P)
  rho.mat[upper.tri(rho.mat, diag = FALSE)] <- rho
  rho.mat <- rho.mat + t(rho.mat)
  diag(rho.mat) <- 1

  if (vec.approx$dist.func.code == "rdist") {
    uvec <- rep(NA, 7)

    uvec[1] <- sig2[ondx[1]]^(-1 / 2)
    uvec[3] <- nugget[ondx[1]]^(-1 / 2)
    uvec[2] <- -1 * uvec[3]
    uvec[7] <- nugget[ondx[2]]^(-1 / 2)
    uvec[6] <- -1 * uvec[7]

    vii <- smoothness[ondx[1]]
    vjj <- smoothness[ondx[2]]
    vij <- (vii + vjj) / 2
    aii <- 1 / rangep[ondx[1]]
    ajj <- 1 / rangep[ondx[2]]
    aij <- sqrt((aii^2 + ajj^2) / 2)
    K1 <- rho.mat[ondx[1], ondx[2]] * sqrt(sig2[ondx[1]]) * sqrt(sig2[ondx[2]]) *
      aii^vii * ajj^vjj * gamma(vij) / (aij^(2 * vij) * sqrt(gamma(vii) *
        gamma(vjj))) *
      cov_func(dist_func(olocs[1, , drop = FALSE], olocs[2, , drop = FALSE], ),
        smoothness = vij, alpha = aij
      )
    K2 <- sig2[ondx[1]]
    bi <- K1 / K2
    ri <- sig2[ondx[2]] - bi * K1

    uvec[5] <- ri^(-1 / 2)
    uvec[4] <- -1 * bi * ri^(-1 / 2)

    vijs <- outer(smoothness, smoothness, "+") / 2
    aijs <- outer(rangep, rangep, function(x, y) sqrt((1 / x^2 + 1 / y^2) / 2))
    gammas <- outer(smoothness, smoothness, function(x, y) gamma((x + y) / 2) / sqrt(gamma(x) * gamma(y)))
    expprod <- outer(rangep^(-smoothness), rangep^(-smoothness), "*")
    sigs <- outer(sig2, sig2, function(x, y) sqrt(x) * sqrt(y))

    full_const <- sigs * gammas * expprod * rho.mat / (aijs^(2 * vijs))

    if (sum(is.na(full_const)) > 0) {
      browser()
    }

    m <- length(q.list$q.y[[n]]) + length(q.list$q.z[[n]])
    q.list$q.y <- lapply(q.list$q.y, function(x) c(x, rep(NA, m - length(x))))
    q.list$q.z <- lapply(q.list$q.z, function(x) c(x, rep(NA, m - length(x))))
    cur.qys <- do.call(cbind, q.list$q.y)
    cur.qzs <- do.call(cbind, q.list$q.z)
    # browser()
    U <- createU_helper_mat(olocs, ondx, cur.qys, cur.qzs, vijs, aijs, full_const, nugget, sig2, uvec)
    # U <- sparseMatrix(i=U1[1,], j=U1[2,], x = U1[3,], triangular = TRUE)
  }

  if (vec.approx$dist.func.code == "custom") {
    # U <- matrix(0, nrow=2*n, ncol=2*n)
    # U <- sparseMatrix(i = 1, j = 1, x = sig2[ondx[1]]^(-1/2), dims = c(2*n, 2*n), triangular = T)
    U1 <- matrix(ncol = 3, nrow = 7)
    U1[1, ] <- c(1, 1, sig2[ondx[1]]^(-1 / 2))
    U1[2, ] <- c(2, 2, nugget[ondx[1]]^(-1 / 2))
    U1[3, ] <- c(1, 2, -U1[2, 3])
    U1[4, ] <- c(4, 4, nugget[ondx[2]]^(-1 / 2))
    U1[5, ] <- c(3, 4, -U1[4, 3])

    # U[1,1] <- sig2[ondx[1]]^(-1/2)
    # U[2,2] <- nugget[ondx[1]]^(-1/2)
    # U[1,2] <- -1*U[2,2]
    # U[4,4] <- nugget[ondx[2]]^(-1/2)
    # U[3,4] <- -1*U[4,4]
    vii <- smoothness[ondx[1]]
    vjj <- smoothness[ondx[2]]
    vij <- (vii + vjj) / 2
    aii <- 1 / rangep[ondx[1]]
    ajj <- 1 / rangep[ondx[2]]
    aij <- sqrt((aii^2 + ajj^2) / 2)
    K1 <- rho.mat[ondx[1], ondx[2]] * sqrt(sig2[ondx[1]]) * sqrt(sig2[ondx[2]]) *
      aii^vii * ajj^vjj * gamma(vij) / (aij^(2 * vij) * sqrt(gamma(vii) *
        gamma(vjj))) *
      cov_func(dist_func(olocs[1, , drop = FALSE], olocs[2, , drop = FALSE], ),
        smoothness = vij, alpha = aij
      )
    K2 <- sig2[ondx[1]]
    bi <- K1 / K2
    ri <- sig2[ondx[2]] - bi * K1
    U1[6, ] <- c(3, 3, ri^(-1 / 2))
    U1[7, ] <- c(1, 3, -1 * bi * ri^(-1 / 2))
    # U[3,3] <- ri^(-1/2)
    # U[1,3] <- -1*bi*ri^(-1/2)
    U2 <- foreach(i = 3:n, .combine = rbind) %dopar% {
      # U[2*i,2*i] <- nugget[ondx[i]]^(-1/2)
      # U[2*i-1,2*i] <- -1*U[2*i,2*i]
      cur.qy <- q.list$q.y[[i]]
      cur.qz <- q.list$q.z[[i]]
      cur.q <- c(cur.qy, cur.qz)
      nq <- length(cur.q)
      cur.U <- matrix(nrow = 3 + nq, ncol = 3)
      # Computing K1 and K2 will be slow because I don't think I can
      # vectorize these calculations unless there is a common smoothness
      # parameter, which generally will not be true in the multivariate
      # case. Rewriting these calculations in C++ is probably going
      # to be the simplest way to speed this up. Currently it is going
      # to be too slow for large data sets.
      K1 <- rep(NA, nq)
      for (j in 1:nq) {
        vii <- smoothness[ondx[i]]
        vjj <- smoothness[ondx[cur.q[j]]]
        vij <- (vii + vjj) / 2
        aii <- 1 / rangep[ondx[i]]
        ajj <- 1 / rangep[ondx[cur.q[j]]]
        aij <- sqrt((aii^2 + ajj^2) / 2)
        # The funky multiplier before the covariance function is necessary
        # is necessary to ensure that the final covariance matrix is
        # positive definite. See equation (9) in Apanasovich (2011).
        K1[j] <- rho.mat[ondx[i], ondx[cur.q[j]]] *
          sqrt(sig2[ondx[i]]) * sqrt(sig2[ondx[cur.q[j]]]) *
          aii^vii * ajj^vjj * gamma(vij) / (aij^(2 * vij) *
            sqrt(gamma(vii) * gamma(vjj))) *
          cov_func(
            dist_func(
              olocs[i, , drop = FALSE],
              olocs[cur.q[j], , drop = FALSE]
            ),
            smoothness = vij, alpha = aij
          )
      }
      K2 <- matrix(nrow = nq, ncol = nq)
      for (j in 1:nq) {
        for (k in j:nq) {
          vii <- smoothness[ondx[cur.q[j]]]
          vjj <- smoothness[ondx[cur.q[k]]]
          vij <- (vii + vjj) / 2
          aii <- 1 / rangep[ondx[cur.q[j]]]
          ajj <- 1 / rangep[ondx[cur.q[k]]]
          aij <- sqrt((aii^2 + ajj^2) / 2)
          K2[j, k] <- rho.mat[ondx[cur.q[j]], ondx[cur.q[k]]] *
            sqrt(sig2[ondx[cur.q[j]]]) * sqrt(sig2[ondx[cur.q[k]]]) *
            aii^vii * ajj^vjj * gamma(vij) /
            (aij^(2 * vij) * sqrt(gamma(vii) * gamma(vjj))) *
            cov_func(
              dist_func(
                olocs[cur.q[j], , drop = FALSE],
                olocs[cur.q[k], , drop = FALSE]
              ),
              smoothness = vij, alpha = aij
            )
          if (j != k) {
            K2[k, j] <- K2[j, k]
          }
        }
      }
      K2 <- K2 + diag(c(rep(0, length(cur.qy)), nugget[ondx[cur.qz]]),
        nrow = nq, ncol = nq
      )
      bi <- t(solve(K2, K1))
      ri <- sig2[ondx[i]] - bi %*% K1
      cur.br <- -1 * as.vector(bi) * c(ri^(-1 / 2))
      cur.U[1, ] <- c(2 * i - 1, 2 * i - 1, ri^(-1 / 2))
      cur.U[2, ] <- c(2 * i, 2 * i, nugget[ondx[i]]^(-1 / 2))
      cur.U[3, ] <- c(2 * i - 1, 2 * i, -cur.U[2, 3])
      temp_len.y <- length(cur.qy)
      # U[(ind+3), ] <- c(2*cur.qy-1,2*i-1, cur.br[1:length(cur.qy)])
      # ind <- ind + 4
      cur.U[4:(4 + temp_len.y - 1), 1] <- 2 * cur.qy - 1
      # show("NO")
      cur.U[4:(4 + temp_len.y - 1), 2] <- 2 * i - 1
      cur.U[4:(4 + temp_len.y - 1), 3] <- cur.br[1:temp_len.y]
      # show("YES")
      # U[2*i-1,2*i-1] <- ri^(-1/2)
      # U[2*cur.qy-1,2*i-1] <- cur.br[1:length(cur.qy)]
      temp_len.z <- length(cur.qz)
      if (temp_len.z > 0) {
        # U[2*cur.qz,2*i-1] <- cur.br[-(1:length(cur.qy))]
        cur.U[(4 + temp_len.y):(3 + nq), 2] <- 2 * i - 1
        cur.U[(4 + temp_len.y):(3 + nq), 1] <- 2 * cur.qz
        cur.U[(4 + temp_len.y):(3 + nq), 3] <- cur.br[-(1:temp_len.y)]
      }
      cur.U
    }
    U <- rbind(U1, U2)
    # browser()
    U <- sparseMatrix(
      i = U[, 1], j = U[, 2], x = U[, 3], dims = c(2 * n, 2 * n),
      triangular = TRUE
    )
  }
  # I think the code below only works for obspred ordering. This probably
  # needs to be fixed when (if) we support other orderings.
  if (sum(!vec.approx$obs) > 0) {
    if (vec.approx$ord.pred != "obspred") {
      stop("Currently only obspred ordering is supported")
    } else {
      drop.seq <- seq(from = (2 * sum(vec.approx$obs) + 2), to = ncol(U), by = 2)
      U <- U[-drop.seq, -drop.seq]
    }
  }
  latent <- rep(TRUE, nrow(U))
  latent[seq(from = 2, to = (2 * sum(vec.approx$obs)), by = 2)] <- FALSE
  return(list(
    U = U, latent = latent, ord = vec.approx$ord, obs = vec.approx$obs,
    ord.pred = vec.approx$ord.pred, ord.z = vec.approx$ord.z,
    cond.yz = vec.approx$cond.yz, ic0 = FALSE
  ))
}
