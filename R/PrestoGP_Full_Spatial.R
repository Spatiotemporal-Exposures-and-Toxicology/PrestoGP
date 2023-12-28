#' Spatial only model created with a likelihood function conditioned on all observations.
#'
#' @slot model SpatialModel
#'
#' @export FullSpatialModel
#'
#' @examples
#' @include PrestoGP_Vecchia_Spatial.R
#' @noRd
FullSpatialModel <- setClass("FullSpatialModel",
  contains = "SpatialModel",
  slots = c(
    model = "SpatialModel"
  ),
  prototype = c(
    model = SpatialModel()
  )
)

validityFullSpatialModel <- function(object) {
  TRUE
}
setValidity("FullSpatialModel", validityFullSpatialModel)
setMethod("initialize", "FullSpatialModel", function(.Object, ...) {
  .Object <- callNextMethod()
  .Object@n_neighbors <- 0
  .Object@min_m <- 0
  validObject(.Object)
  .Object
})

#' specify
#'
#' Do nothing, full models do not specify a conditioning set, they condition on all observations
#'
#' @param model A model object.
#' @param locs a locations matrix
#' @param m the number of neighbors to condition on (will always be n-1 for full models)
#'
#' @return the unmodified model
setMethod("specify", "FullSpatialModel", function(model, locs, m) {
  invisible(model)
})

setMethod("compute_residuals", "FullSpatialModel", function(model, Y, Y.hat) {
  model@res <- as.double(Y - Y.hat)
  invisible(model)
})

setMethod("estimate_theta", "FullSpatialModel", function(model, locs) {
  n <- length(model@Y_train)
  full.result <- optim(
    par = log(model@covparams), fn = negloglik_full_spatial,
    y = model@res, locs = locs, N = n, method = "Nelder-Mead",
    control = list(trace = 0)
  )
  model@covparams <- exp(full.result$par)
  model@LL_Vecchia_krig <- full.result$value
  invisible(model)
})

#' transform_data
#'
#' Model transformation for full models is performed directly on X and Y matrices using the inverted Cholesky decompostion of the covariance matrix.
#' (no Vecchia approximation)
#'
#' @param model A model object.
#' @param Y the dependent variable matrix
#' @param X the independent variable matrix
#'
#' @return the unmodified model
setMethod("transform_data", "FullSpatialModel", function(model, Y, X) {
  n <- length(model@Y_train)
  Sigma.oo <- model@covparams[1] * Exponential(rdist(model@locs_train), range = model@covparams[2]) + model@covparams[3] * diag(n)
  Omega.lc <- solve(t(chol(Sigma.oo)))
  model@y_tilde <- Matrix(Omega.lc %*% Y)
  model@X_tilde <- Matrix(Omega.lc %*% X)
  invisible(model)
})

setMethod("theta_names", "FullSpatialModel", function(model) {
  c("Marginal Variance", "Spatial Range", "Nugget")
})
