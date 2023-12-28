#' Spatiotemporal model created with a likelihood function conditioned on a subset of the observations.
#'
#' @slot model PrestoGPModel.
#'
#' @export SpatiotemporalModel
#'
#' @examples
#' @include PrestoGP_Model.R
#' @noRd
SpatiotemporalModel <- setClass("SpatiotemporalModel",
  contains = "PrestoGPModel",
  slots = c(
    model = "PrestoGPModel"
  )
)

validitySpatiotemporalModel <- function(object) {
  TRUE
}
setValidity("SpatiotemporalModel", validitySpatiotemporalModel)

setMethod("initialize", "SpatiotemporalModel", function(.Object, n_neighbors = 25, ...) {
  .Object@n_neighbors <- n_neighbors
  .Object@min_m <- 3
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
setMethod("prestogp_predict", "SpatiotemporalModel", function(model, X, locs, m = NULL) {
  # validate parameters
  if (!is.matrix(X)) {
    stop("X parameter must be a matrix.")
  }
  if (!is.matrix(locs)) {
    stop("The locs parameter must be a matrix.")
  }
  if (nrow(X) != nrow(locs)) {
    stop("The number of locations must match the number of X observations.")
  }
  if (ncol(locs) != 3) {
    stop("The locs parameter must have 3 columns.")
  }
  if (is.null(m)) { # m defaults to the value used for training
    m <- model@n_neighbors
  }
  if (m == 0) {
    m <- nrow(model@X_train) - 1
  }

  # Vecchia prediction at new locations
  # Vecchia.Pred <- predict(model@Vecchia_SCAD_fit[[1]], X = X, which = model@lambda_1se_idx[[1]])
  Vecchia.Pred <- predict(model@linear_model, newx = X, s = model@linear_model$lambda[model@lambda_1se_idx])
  # Vecchia trend prediction at observed data
  # Vecchia.hat <- predict(model@Vecchia_SCAD_fit[[1]], X = model@X_train, which = model@lambda_1se_idx[[1]])
  Vecchia.hat <- predict(model@linear_model, newx = model@X_train, s = model@linear_model$lambda[model@lambda_1se_idx])

  # Test set prediction
  res <- model@Y_train - Vecchia.hat

  locs.train.scaled <- scale_locs(model, model@locs_train)
  locs.scaled <- scale_locs(model, locs)
  vec.approx.test <- vecchia_specify(locs.train.scaled, m, locs.pred = locs.scaled)

  ## carry out prediction
  pred <- vecchia_prediction(res, vec.approx.test, c(model@covparams[1], 1, 0.5), model@covparams[4])

  # prediction function can return both mean and sds
  # returns a list with elements mu.pred,mu.obs,var.pred,var.obs,V.ord
  Vec.mean <- pred$mu.pred + Vecchia.Pred # residual + mean trend
  # option to include or exclude theta below
  Vec.sds <- sqrt(pred$var.pred + model@covparams[4]) # standard deviation

  return(list("means" = Vec.mean, "standard deviations" = Vec.sds))
})

#' calc_covparams
#'
#' Set initial value of covarariance parameters.
#'
#' @param model The model to set the covariance parameters of
#' @param locs the locations matrix
#' @param Y the dependent variable matrix
#'
#' @return a model with initial covariance parameters
setMethod("calc_covparams", "SpatiotemporalModel", function(model, locs, Y) {
  N <- length(Y)
  d.sample <- sample(1:N, max(2, ceiling(N / 50)), replace = FALSE) # sample(1:N,N/50,replace = FALSE)
  D.sample <- rdist(locs[d.sample, 1:2])
  t.sample <- rdist(locs[d.sample, 3])
  model@covparams <- c(.9 * var(Y), mean(D.sample) / 4, mean(t.sample) / 4, 0.1 * var(Y)) # smoothness=0.5
  invisible(model)
})

#' specify
#'
#' Specify the conditioning set using m nearest neighbors.
#'
#' @param model The model to specify
#' @param locs the locations matrix
#' @param m the number of neighbors to condition on
#'
#' @return a model with a specified conditioning set
setMethod("specify", "SpatiotemporalModel", function(model, locs, m) {
  locs.scaled <- scale_locs(model, locs)
  model@vecchia_approx <- vecchia_specify(locs.scaled, m)
  invisible(model)
})

#' scale_locs
#'
#' Scale the locations matrix by the covariance parameters
#'
#' @param model The model with locations to scale
#' @param locs the locations matrix
#'
#' @return a matrix with scaled locations
setMethod("scale_locs", "SpatiotemporalModel", function(model, locs) {
  cbind(locs[, 1] / model@covparams[2], locs[, 2] / model@covparams[2], locs[, 3] / model@covparams[3])
})

#' compute_residuals
#'
#' Compute residuals based on beta parameters
#'
#' @param model The model to compute the residual of
#' @param Y the training dependent variable matrix
#' @param Y.hat the predicted training dependent variable matrix based on beta hat
#'
#' @return a model with computed residuals
setMethod("compute_residuals", "SpatiotemporalModel", function(model, Y, Y.hat) {
  model@res <- as.double(Y - Y.hat)
  model@vecchia_approx$zord <- model@res[model@vecchia_approx$ord]
  invisible(model)
})


#' estimate_theta
#'
#' Estimate covariance parameters during a single round of model optimization
#'
#' @param model The model to estimate theta for
#' @param locs the locations matrix
#'
#' @return a model with an updated covariance parameters estimate
setMethod("estimate_theta", "SpatiotemporalModel", function(model, locs) {
  ord_locs <- locs[model@vecchia_approx$ord, ]
  locs_mat <- cbind(ord_locs[, 1], ord_locs[, 2], ord_locs[, 3])
  vecchia.result <- optim(
    par = log(model@covparams), fn = negloglik_vecchia_ST,
    locs = locs_mat, res = model@res, vecchia.approx = model@vecchia_approx, method = "Nelder-Mead",
    control = list(trace = 0)
  )
  model@LL_Vecchia_krig <- vecchia.result$value
  model@covparams <- exp(vecchia.result$par)
  invisible(model)
})

#' transform_data
#'
#' Transform data to be independent, identically distributed
#'
#' @param model The model to estimate theta for
#' @param Y the dependent variable matrix
#' @param X the independent variable matrix
#'
#' @return a model with i.i.d. data
setMethod("transform_data", "SpatiotemporalModel", function(model, Y, X) {
  transformed.data <- transform_iid(cbind(Y, as.matrix(X)),
    vecchia.approx = model@vecchia_approx, c(model@covparams[1], 1, 0.5), model@covparams[4]
  )
  model@y_tilde <- Matrix(transformed.data[, 1])
  model@X_tilde <- Matrix(transformed.data[, -1])
  invisible(model)
})

#' theta_names
#'
#' Return a vector specifying the names of the different covariance parameters (used for show method)
#'
#' @param model The model to estimate theta for
#'
#' @return a vector with the namems of the covariance parameterss
setMethod("theta_names", "SpatiotemporalModel", function(model) {
  c("Marginal Variance", "Spatial Range", "Temporal Range", "Nugget")
})
