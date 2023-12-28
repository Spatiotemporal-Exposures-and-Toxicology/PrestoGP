#' Spatiotemporal model created with a likelihood function conditioned on all observations.
#'
#' @slot model SpatiotemporalModel
#'
#' @export SpatiotemporalFullModel
#'
#' @examples
#' @include PrestoGP_Vecchia_Spatiotemporal.R
#' @noRd
SpatiotemporalFullModel <- setClass("SpatiotemporalFullModel",
  contains = "SpatiotemporalModel",
  slots = c(
    model = "SpatiotemporalModel"
  )
)

validitySpatiotemporalFullModel <- function(object) {
  TRUE
}
setValidity("SpatiotemporalFullModel", validitySpatiotemporalFullModel)
setMethod("initialize", "SpatiotemporalFullModel", function(.Object, ...) {
  .Object <- callNextMethod()
  .Object@n_neighbors <- 0
  .Object@min_m <- 0
  validObject(.Object)
  .Object
})

setMethod("calc_covparams", "SpatiotemporalFullModel", function(model, locs, Y) {
  N <- length(Y)
  d.sample <- sample(1:N, max(2, ceiling(N / 50)), replace = FALSE)
  D.sample <- rdist(locs[d.sample, 1:2])
  t.sample <- rdist(locs[d.sample, 3])
  model@covparams <- c(.9 * var(Y), mean(D.sample) / 4, mean(t.sample) / 4, 0.1 * var(Y))
  invisible(model)
})

setMethod("specify", "SpatiotemporalFullModel", function(model, locs, m) {
  invisible(model)
})

setMethod("compute_residuals", "SpatiotemporalFullModel", function(model, Y, Y.hat) {
  model@res <- as.double(Y - Y.hat)
  # model@vecchia_approx$zord = model@res[model@vecchia_approx$ord]
  invisible(model)
})

setMethod("estimate_theta", "SpatiotemporalFullModel", function(model, locs) {
  n <- length(model@Y_train)
  full.result <- optim(
    par = log(model@covparams), fn = negloglik_full_ST,
    y = model@res, locs = locs, N = n, method = "Nelder-Mead",
    control = list(trace = 0)
  )
  model@covparams <- exp(full.result$par)
  model@LL_Vecchia_krig <- full.result$value
  invisible(model)
})

setMethod("transform_data", "SpatiotemporalFullModel", function(model, Y, X) {
  n <- length(model@Y_train)
  locs.scaled <- scale_locs(model, model@locs_train)
  Omega.full <- model@covparams[1] * Exponential(rdist(locs.scaled), range = 1) + model@covparams[4] * diag(n)
  Omega.lc <- solve(t(chol(Omega.full)))
  model@y_tilde <- Matrix(Omega.lc %*% Y)
  model@X_tilde <- Matrix(Omega.lc %*% X)
  invisible(model)
})
