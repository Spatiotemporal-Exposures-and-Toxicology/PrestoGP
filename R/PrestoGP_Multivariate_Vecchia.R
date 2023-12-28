#' Title
#'
#' @slot model PrestoGPModel.
#'
#' @export MultivariateVecchiaModel
#'
#' @examples
#' @include PrestoGP_Model.R
#' @noRd
MultivariateVecchiaModel <- setClass("MultivariateVecchiaModel",
  contains = "PrestoGPModel",
  slots = c(
    model = "PrestoGPModel",
    param_sequence = "matrix",
    logparams = "numeric"
  )
)

validityMultivariateVecchiaModel <- function(object) {
  TRUE
}
setValidity(
  "MultivariateVecchiaModel",
  validityMultivariateVecchiaModel
)

setMethod("initialize", "MultivariateVecchiaModel", function(.Object, ...) {
  .Object@n_neighbors <- 0 # 25
  .Object@min_m <- 0
  .Object <- callNextMethod()
  validObject(.Object)
  .Object
})

#' Make predictions using a previously fit spatiotemporal model.
#'
#' @param model A model object returned by prestogp_fit.
#' @param X Independent variable matrix for prediction.
#' @param locs Locations matrix for which the value is predicted.
#' @param m The number of neighbors to condition on (if not provided, training value of m will be used).
#'
#' @return predicted values for the dependent variable
#' @export
#'
#' @examples
#'
#' ...
#' model <- SpatiotemporalModel()
#' model <- prestogp_fit(model, logNO2, X, locs)
#' prediction <- prestogp_predict(model, X.test, locs.test)
#' Vec.mean <- prediction[[1]]
#' Vec.sds <- prediction[[2]]
setMethod("prestogp_predict", "MultivariateVecchiaModel", function(model, X, locs, m = NULL, ordering.pred = c("obspred", "general"), pred.cond = c("independent", "general"), return.values = c("mean", "meanvar")) {
  # validate parameters
  ordering.pred <- match.arg(ordering.pred)
  pred.cond <- match.arg(pred.cond)
  return.values <- match.arg(return.values)
  if (!is.list(X)) {
    stop("X parameter must be a list.")
  }
  if (!is.list(locs)) {
    stop("locs parameter must be a list.")
  }
  ndx.out <- NULL
  for (i in 1:length(locs)) {
    if (nrow(X[[i]]) != nrow(locs[[i]])) {
      stop("The number of locations must match the number of X observations.")
    }
    ndx.out <- c(ndx.out, rep(i, nrow(locs[[i]])))
  }
  X <- psych::superMatrix(X)
  if (ncol(X) != ncol(model@X_train)) {
    stop("The number of predictors in X must match the training data")
  }
  if (is.null(m)) { # m defaults to the value used for training
    m <- model@n_neighbors
  }
  stopifnot((m > 0)) # FIXME m is not required by full model

  # Vecchia prediction at new locations
  Vecchia.Pred <- predict(model@linear_model, newx = X, s = model@linear_model$lambda[model@lambda_1se_idx])
  # Vecchia trend prediction at observed data
  Vecchia.hat <- predict(model@linear_model, newx = model@X_train, s = model@linear_model$lambda[model@lambda_1se_idx])

  # Test set prediction
  res <- model@Y_train - Vecchia.hat

  locs.train.scaled <- scale_locs(model, model@locs_train)
  locs.scaled <- scale_locs(model, locs)
  vec.approx.test <- vecchia_Mspecify(locs.train.scaled, m,
    locs.list.pred = locs.scaled,
    ordering.pred = ordering.pred,
    pred.cond = pred.cond
  )

  ## carry out prediction
  if (!model@apanasovich) {
    params <- model@covparams
    param.seq <- model@param_sequence
    pred <- vecchia_Mprediction(res, vec.approx.test,
      c(
        params[1:param.seq[1, 2]],
        rep(1, param.seq[2, 2] - param.seq[2, 1] + 1),
        params[param.seq[3, 1]:
        param.seq[5, 2]]
      ),
      return.values = return.values
    )
  } else {
    pred <- vecchia_Mprediction(res, vec.approx.test, model@covparams,
      return.values = return.values
    )
  }

  # prediction function can return both mean and sds
  # returns a list with elements mu.pred,mu.obs,var.pred,var.obs,V.ord
  Vec.mean <- pred$mu.pred + Vecchia.Pred # residual + mean trend
  if (return.values == "mean") {
    return.list <- list(means = Vec.mean)
  } else {
    warning("Variance estimates do not include model fitting variance and are anticonservative. Use with caution.")
    vec.sds <- sqrt(pred$var.pred)
    for (i in 1:length(locs)) {
      vec.sds[ndx.out == i] <- sqrt(vec.sds[ndx.out == i] +
        model@covparams[model@param_sequence[4, i]])
    }
    return.list$sds <- vec.sds
  }
  # option to include or exclude theta below
  #  Vec.sds = sqrt(pred$var.pred + model@covparams[4]) #standard deviation

  return(return.list)
})

setMethod("calc_covparams", "MultivariateVecchiaModel", function(model, locs, Y) {
  # cor.matrix <- cor(Y)
  # P <- ncol(Y)
  # col.vars <- apply(Y, 2, var)
  # N <- length(Y)
  P <- length(Y)
  col.vars <- rep(NA, P)
  D.sample.bar <- rep(NA, model@nscale * P)
  for (i in 1:P) {
    col.vars[i] <- var(Y[[i]])
    N <- length(Y[[i]])
    # TODO find a better way to compute initial spatial range
    for (j in 1:model@nscale) {
      d.sample <- sample(1:N, max(2, ceiling(N / 50)), replace = FALSE)
      D.sample <- rdist(locs[[i]][d.sample, model@scaling == j])
      D.sample.bar[(i - 1) * model@nscale + j] <- mean(D.sample) / 4
    }
  }
  model@logparams <- create.initial.values.flex(
    c(0.9 * col.vars), # marginal variance
    D.sample.bar, # range
    rep(0.5, P), # smoothness
    c(.1 * col.vars), # nuggets
    rep(0, choose(P, 2)),
    P
  )
  model@param_sequence <- create.param.sequence(P, model@nscale)
  model <- transform_covariance_parameters(model)
  invisible(model)
})

setMethod("specify", "MultivariateVecchiaModel", function(model, locs, m) {
  locs.scaled <- scale_locs(model, locs)
  model@vecchia_approx <- vecchia_Mspecify(locs.scaled, m)
  if (!model@apanasovich) {
    olocs.scaled <- model@vecchia_approx$locsord
    for (i in 1:length(locs)) {
      for (j in 1:model@nscale) {
        olocs.scaled[model@vecchia_approx$ondx == i, model@scaling == j] <-
          olocs.scaled[model@vecchia_approx$ondx == i, model@scaling == j] *
            model@covparams[model@param_sequence[2, 1] +
              model@nscale * (i - 1) + j - 1]
      }
    }
    model@vecchia_approx$locsord <- olocs.scaled
  }
  invisible(model)
})

setMethod("scale_locs", "MultivariateVecchiaModel", function(model, locs) {
  if (model@apanasovich) {
    return(locs)
  } else {
    locs.out <- locs
    for (i in 1:length(locs)) {
      for (j in 1:model@nscale) {
        locs.out[[i]][, model@scaling == j] <-
          locs[[i]][, model@scaling == j] /
            model@covparams[model@param_sequence[2, 1] +
              model@nscale * (i - 1) + j - 1]
      }
    }
    return(locs.out)
  }
})

setMethod("compute_residuals", "MultivariateVecchiaModel", function(model, Y, Y.hat) {
  model@res <- as.double(Y - Y.hat)
  model@vecchia_approx$zord <- model@res[model@vecchia_approx$ord]
  invisible(model)
})

setMethod("estimate_theta", "MultivariateVecchiaModel", function(model, locs, optim.control, method) {
  P <- length(locs)
  #  locs_list <- list()
  #  y <- list()
  #  for(i in 1:P){
  #    locs_list[[i]] <- locs
  #    y[[i]] <- model@Y_train[,i]
  #  }
  #  show(model@covparams)
  #  show(model@param_sequence)
  #  show(model@logparams)
  if (model@apanasovich) {
    vecchia.result <- optim(
      par = model@logparams,
      fn = mvnegloglik,
      vecchia.approx = model@vecchia_approx,
      y = model@res,
      P = P,
      param.seq = model@param_sequence,
      method = method,
      control = optim.control
    )
  } else {
    vecchia.result <- optim(
      par = model@logparams,
      fn = mvnegloglik_ST,
      vecchia.approx = model@vecchia_approx,
      y = model@res,
      P = P,
      param.seq = model@param_sequence,
      scaling = model@scaling,
      nscale = model@nscale,
      method = method,
      control = optim.control
    )
  }

  model@LL_Vecchia_krig <- vecchia.result$value
  model@logparams <- vecchia.result$par
  model <- transform_covariance_parameters(model)
  invisible(model)
})

setMethod("transform_covariance_parameters", "MultivariateVecchiaModel", function(model) {
  P <- length(model@Y_train)
  if (P > 1) {
    model@covparams <- c(
      exp(model@logparams[1:model@param_sequence[2, 2]]),
      gtools::inv.logit(
        model@logparams[model@param_sequence[3, 1]:
        model@param_sequence[3, 2]],
        0, 2.5
      ),
      exp(model@logparams[model@param_sequence[4, 1]:
      model@param_sequence[4, 2]]),
      tanh(model@logparams[model@param_sequence[5, 1]:
      model@param_sequence[5, 2]])
    )
  } else {
    model@covparams <- c(
      exp(model@logparams[1:model@param_sequence[2, 2]]),
      gtools::inv.logit(
        model@logparams[model@param_sequence[3, 1]:
        model@param_sequence[3, 2]],
        0, 2.5
      ),
      exp(model@logparams[model@param_sequence[4, 1]:
      model@param_sequence[4, 2]]), 1
    )
  }
  invisible(model)
})

setMethod("transform_data", "MultivariateVecchiaModel", function(model, Y, X) {
  vecchia.approx <- model@vecchia_approx
  if (!model@apanasovich) {
    params <- model@covparams
    param.seq <- model@param_sequence
    olocs.scaled <- vecchia.approx$locsord
    for (i in 1:vecchia.approx$P) {
      for (j in 1:model@nscale) {
        olocs.scaled[model@vecchia_approx$ondx == i, model@scaling == j] <-
          olocs.scaled[
            model@vecchia_approx$ondx == i,
            model@scaling == j
          ] /
            model@covparams[param.seq[2, 1] + model@nscale * (i - 1) + j - 1]
      }
    }
    vecchia.approx$locsord <- olocs.scaled
    transformed.data <- transform_miid(cbind(Y, as.matrix(X)),
      vecchia.approx = vecchia.approx,
      c(
        params[1:param.seq[1, 2]],
        rep(
          1,
          param.seq[2, 2] - param.seq[2, 1] + 1
        ),
        params[param.seq[3, 1]:
        param.seq[5, 2]]
      )
    )
  } else {
    transformed.data <- transform_miid(cbind(Y, as.matrix(X)),
      vecchia.approx = vecchia.approx,
      model@covparams
    )
  }

  xcols <- ncol(model@X_train)
  ycols <- ncol(model@Y_train)
  tcols <- ncol(transformed.data)
  model@y_tilde <- Matrix(transformed.data[, 1:ncol(model@Y_train)])
  model@X_tilde <- Matrix(transformed.data[, (ncol(model@Y_train) + 1):ncol(transformed.data)], sparse = FALSE)
  invisible(model)
})

setMethod("theta_names", "MultivariateVecchiaModel", function(model) {
  c("Marginal Variance", "Range", "Smoothness", "Nugget")
})
