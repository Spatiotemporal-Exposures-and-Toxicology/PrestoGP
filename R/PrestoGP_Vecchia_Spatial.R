#' Spatial only model created with a likelihood function conditioned on a subset of the observations.
#'
#' @slot model PrestoGPModel
#'
#' @export SpatialModel
#'
#' @examples
#' @include PrestoGP_Model.R
#' @noRd
SpatialModel <- setClass("SpatialModel",
  contains = "PrestoGPModel",
  slots = c(
    model = "PrestoGPModel"
  )
)

validitySpatialModel <- function(object) {
  if (object@n_neighbors < object@min_m) {
    stop(paste("N_neighbors must be at least ", object@min_m, ".", sep = ""))
  }
  TRUE
}
setValidity("SpatialModel", validitySpatialModel)

setMethod("initialize", "SpatialModel", function(.Object, n_neighbors = 25, ...) {
  .Object@n_neighbors <- n_neighbors
  .Object@min_m <- 3
  .Object <- callNextMethod()
  # validObject(.Object)
  .Object
})

#' Make predictions using a previously fit spatial model.
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
setMethod("prestogp_predict", "SpatialModel", function(model, X, locs, m = NULL) {
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
  if (ncol(locs) != 2) {
    stop("The locs parameter must have 2 columns.")
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
  pred <- vecchia_prediction(res, vec.approx.test, c(model@covparams[1], model@covparams[2], 0.5), model@covparams[3])
  # prediction function can return both mean and sds
  # returns a list with elements mu.pred,mu.obs,var.pred,var.obs,V.ord
  Vec.mean <- pred$mu.pred + Vecchia.Pred # residual + mean trend
  # option to include or exclude theta below
  # Vec.sds = sqrt(pred$var.pred + model@covparams[4]) #standard deviation
  Vec.sds <- sqrt(pred$var.pred)
  return(list("means" = Vec.mean, "standard deviations" = Vec.sds))
})

setMethod("calc_covparams", "SpatialModel", function(model, locs, Y) {
  N <- length(Y)
  d.sample <- sample(1:N, max(2, ceiling(N / 50)), replace = FALSE)
  D.sample <- rdist(locs[d.sample, 1:2])
  model@covparams <- c(.9 * var(Y), mean(D.sample) / 4, 0.1 * var(Y))
  invisible(model)
})

setMethod("specify", "SpatialModel", function(model, locs, m) {
  model@vecchia_approx <- vecchia_specify(locs, m)
  invisible(model)
})

setMethod("scale_locs", "SpatialModel", function(model, locs) {
  locs
})

setMethod("compute_residuals", "SpatialModel", function(model, Y, Y.hat) {
  model@res <- as.double(Y - Y.hat)
  model@vecchia_approx$zord <- model@res[model@vecchia_approx$ord]
  invisible(model)
})

setMethod("estimate_theta", "SpatialModel", function(model, locs, optim.control) {
  vecchia.result <- optim(
    par = log(model@covparams), fn = negloglik_vecchia,
    locs = locs, res = model@res, vecchia.approx = model@vecchia_approx, method = "Nelder-Mead",
    control = optim.control
  )
  model@LL_Vecchia_krig <- vecchia.result$value
  model@covparams <- exp(vecchia.result$par)
  invisible(model)
})

setMethod("transform_data", "SpatialModel", function(model, Y, X) {
  transformed.data <- transform_iid(cbind(Y, as.matrix(X)),
    vecchia.approx = model@vecchia_approx, c(model@covparams[1], model@covparams[2], 0.5), model@covparams[3]
  )
  model@y_tilde <- Matrix(transformed.data[, 1])
  model@X_tilde <- Matrix(transformed.data[, -1])
  invisible(model)
})

setMethod("theta_names", "SpatialModel", function(model) {
  c("Marginal Variance", "Spatial Range", "Nugget")
})
