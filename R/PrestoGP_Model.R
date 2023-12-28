# Look at https://github.com/cran/gpclib/blob/master/R/Rgpc.R
library(methods)

setOldClass("cv.glmnet")


#' PrestoGPModel
#'
#' @slot covparams numeric.
#' @slot beta numeric.
#' @slot lambda_1se_idx numeric.
#' @slot vecchia_approx list.
#' @slot y_tilde numeric.
#' @slot X_tilde dgeMatrix.
#' @slot res numeric.
#' @slot Vecchia_SCAD_fit cv.ncvreg.
#'
#' @export PrestoGPModel
#' @import GPvecchia Matrix fields ncvreg readxl scoringRules MASS glmnet
#' @importFrom stats optim predict var
#' @importFrom aod wald.test
#' @importFrom dplyr %>%
#'
#' @examples
#' @noRd
PrestoGPModel <- setClass(
  "PrestoGPModel",
  slots = list(
    covparams = "numeric",
    beta = "matrix",
    # "numeric", #the beta matrix
    lambda_1se_idx = "numeric",
    # "list", #the index of the best model
    vecchia_approx = "list",
    # the output of vecchia_specify
    y_tilde = "dgeMatrix",
    # iid transformed dependent variables matrix
    X_tilde = "dgeMatrix",
    # iid transformed independent variables matrix
    res = "numeric",
    # residuals
    linear_model = "cv.glmnet",
    # the linear model
    X_train = "matrix",
    # the original independent variable matrix
    Y_train = "matrix",
    # the original dependent variable matrix
    locs_train = "list",
    # the location / temporal matrix
    converged = "logical",
    # a logical variable that is true if
    # the model fitting process has converged
    LL_Vecchia_krig = "numeric",
    # the value of the negative log
    # likelihood function after optimization
    error = "numeric",
    # negative log likelihood + SCAD penalty likelihood
    n_neighbors = "numeric",
    # the number of neighbors to condition on for
    # the Vecchia approximation
    min_m = "numeric",
    # the minimum m required by the specific model type
    # (full vs Vecchia)
    alpha = "numeric",
    # the alpha ratio of ridge to lasso penalty
    scaling = "numeric",
    # the indices of the scale parameters,
    nscale = "numeric",
    # the number of scale parameters
    apanasovich = "logical"
  ) # should the Apanasovich model be used
)

validityPrestoGPModel <- function(object) {
  TRUE
}
setValidity("PrestoGPModel", validityPrestoGPModel)

setMethod("initialize", "PrestoGPModel", function(.Object, ...) {
  .Object@linear_model <- structure(list(), class = "cv.glmnet")
  .Object@alpha <- 1 # 0.5
  .Object <- callNextMethod()
  validObject(.Object)
  .Object
})

setGeneric("show_theta", function(object, Y_names) {
  standardGeneric("show_theta")
})
setGeneric("prestogp_fit", function(model,
                                    Y,
                                    X,
                                    locs,
                                    scaling = NULL,
                                    apanasovich = FALSE,
                                    covparams = NULL,
                                    beta.hat = NULL,
                                    tol = 0.999999,
                                    max_iters = 100,
                                    verbose = FALSE,
                                    optim.method = "Nelder-Mead",
                                    optim.control = list(
                                      trace = 0,
                                      reltol = 1e-4,
                                      maxit = 5000
                                    ),
                                    parallel = FALSE) {
  standardGeneric("prestogp_fit")
})
setGeneric("prestogp_predict", function(model,
                                        X = "matrix",
                                        locs = "matrix",
                                        m = "numeric",
                                        ordering.pred = c("obspred", "general"),
                                        pred.cond = c("independent", "general"),
                                        return.values = c("mean", "meanvar")) {
  standardGeneric("prestogp_predict")
})
setGeneric("calc_covparams", function(model, locs, Y) {
  standardGeneric("calc_covparams")
})
setGeneric("specify", function(model, locs, m) {
  standardGeneric("specify")
})
setGeneric("compute_residuals", function(model, Y, Y.hat) {
  standardGeneric("compute_residuals")
})
setGeneric("transform_data", function(model, Y, X) {
  standardGeneric("transform_data")
})
setGeneric("estimate_theta", function(model, locs, optim.control, method) {
  standardGeneric("estimate_theta")
})
setGeneric("estimate_betas", function(model, parallel) {
  standardGeneric("estimate_betas")
})
setGeneric("compute_error", function(model, y, X) {
  standardGeneric("compute_error")
})
setGeneric("scale_locs", function(model, locs) {
  standardGeneric("scale_locs")
})
setGeneric("theta_names", function(model) {
  standardGeneric("theta_names")
})
setGeneric("transform_covariance_parameters", function(model) {
  standardGeneric("transform_covariance_parameters")
})


#' show
#'
#' Print a summary of the model and its parameters
#'
#' @param object the PrestoGP model object
setMethod(
  "show", "PrestoGPModel", # TODO consider exporting this
  function(object) {
    cat("PrestoGP Model\n") # TODO print out type of model
    cat("Negative Log-Likelihood: ", object@error, "\n")
    cat("Covariance Parameters:\n")
    Y_names <- colnames(object@Y_train)
    if (is.null(Y_names)) {
      Y_names <- unlist(lapply(
        seq_len(ncol(object@Y_train)),
        function(x) {
          paste("Outcome", x)
        }
      ))
    }
    show_theta(object, Y_names)

    y_hat <-
      matrix(
        predict(object@linear_model, newx = object@X_train),
        nrow = nrow(object@X_train),
        ncol = ncol(object@Y_train)
      )
    mse <-
      crossprod((object@Y_train - y_hat)) / (nrow(object@Y_train) -
        colSums(object@beta))
    cat("\nTraining MSE: ", diag(mse), "\n")
    X <- cbind(1, object@X_train)
    covm <- MASS::ginv(t(X) %*% X)
    cat("Non-zero Beta Parameters:\n")
    # TODO compare to zero within a tolerance
    # nnz_betas <- lapply(object@beta, 2, function(x){which(x != 0.0)})
    nnz_betas <- list()
    for (col in seq_len(ncol(object@Y_train))) {
      nnz_betas <-
        append(nnz_betas, list(which(object@beta[, col] != 0.0)))
    }
    X_names <- colnames(object@X_train)
    if (is.null(X_names)) {
      X_names <- unlist(lapply(
        seq_len(ncol(object@X_train)),
        function(x) {
          paste("Ind. Variable", x)
        }
      ))
    }
    X_names <- append("Intercept", X_names)
    for (i in seq_len(ncol(object@Y_train))) {
      cat(Y_names[i], " Parameters:\n")
      beta_summary <-
        data.frame(matrix(
          ncol = 4,
          nrow = 0,
          dimnames = list(
            NULL,
            c(
              "Parameter", "Estimate",
              "Standard Error", "Walds P-value"
            )
          )
        ))
      # for(nnz in nnz_betas[i]){
      for (j in seq_along(nnz_betas[[i]])) {
        nnz <- nnz_betas[[i]][[j]]
        walds <-
          wald.test(covm * mse[i, i], object@beta[, i], Terms = nnz)
        std_err <- sqrt(diag(covm) * mse[i, i])
        walds_p <- walds$result$chi2[3]
        beta_summary[nrow(beta_summary) + 1, ] <-
          list(X_names[nnz], object@beta[nnz, i], std_err[nnz], walds_p)
      }
      print(beta_summary, row.names = FALSE)
      cat("\n")
    }
    invisible(object)
  }
)

#' show_theta
#'
#' Print the covariance parameters in a table
#'
#' @param object the PrestoGP model object
#' @param Y_names the names of the different outcome variables
#' (may just be numbers if not provided in training input)
setMethod(
  "show_theta", "PrestoGPModel",
  function(object, Y_names) {
    theta_name_arr <- theta_names(object)
    theta_summary <-
      data.frame(matrix(
        ncol = ncol(object@Y_train) + 1,
        nrow = length(theta_name_arr),
        dimnames = list(NULL, c("Parameter", Y_names))
      ))
    for (i in seq_along(theta_name_arr)) {
      theta_row <-
        object@covparams[((i - 1) * ncol(object@Y_train) + 1):
        (i * ncol(object@Y_train))]
      for (j in seq_len(ncol(object@Y_train))) {
        theta_summary[i, j + 1] <- theta_row[j]
      }
    }
    for (j in seq_along(theta_name_arr)) {
      theta_summary[j, 1] <- theta_name_arr[j]
    }
    print(theta_summary, row.names = FALSE)
    # TODO show Rho matrix if there are 2 or more outcomes
  }
)

#' Train a PrestoGP model.
#'
#' This method fits any PrestoGP model given a matrix of locations, a matrix of
#' #' independent variables, and a matrix of dependent variables.
#'
#' @param model The model object being fit.
#' @param Y A matrix containing training values for the dependent variable.
#' @param X A matrix containing training values for the independent varaibles.
#' @param locs A matrix containing the training spatial coordinates and times.
#' @param covparams The initial covariance parameters to use (optional).
#' @param beta.hat The initial beta parameters to use (optional).
#' @param tol The model is considered converged when error isn't less than
#' tol*previous_error (optional).
#' @param m The number of neighboring datapoints to condition on in the
#' likelihood function (optional).
#' @param verbose If TRUE, additional information about model fit
#'  will be printed.
#'
#' @return An object containing model parameters for spatiotemporal prediction.
#' @export
#'
#' @examples
#'
#' US_NO2_Data <- read.csv("data/US_ST_NO2_Data.csv")
#' lat <- US_NO2_Data$Latitude # Latitude
#' lon <- US_NO2_Data$Longitude # Longitude
#' logNO2 <- log(US_NO2_Data$Y) # ozone data
#' logNO2 <- logNO2[1:1000, ]
#' time <- US_NO2_Data$YearFrac # time in fractional years
#' N <- length(logNO2)
#' locs <- cbind(lon, lat, time) # Coordinates in R3 (x,y,t)
#' locs <- locs[1:1000, ]
#' X.Design <- US_NO2_Data[, 7:145] # Design matrix
#' X.Design <- scale(X.Design)
#' X <- X.Design[1:1000, ]
#' model <- PrestoGPSpatiotemporalModel()
#' model <- prestogp_fit(model, logNO2, X, locs)
#' ...
setMethod(
  "prestogp_fit", "PrestoGPModel",
  function(model,
           Y,
           X,
           locs,
           scaling = NULL,
           apanasovich = FALSE,
           covparams = NULL,
           beta.hat = NULL,
           tol = 0.999999,
           max_iters = 100,
           verbose = FALSE,
           optim.method = "Nelder-Mead",
           optim.control = list(
             trace = 0,
             reltol = 1e-4,
             maxit = 5000
           ),
           parallel = FALSE) {
    # parameter validation
    # TODO: This method should check for input errors in the
    # multivariate case (where Y, X, and locs are lists)
    if (!is.matrix(locs) && !is.list(locs)) {
      stop("locs parameter must be a matrix or a list.")
    }
    if (is.double(Y) && length(Y) == nrow(locs)) {
      Y <- as.matrix(Y)
    }
    if (!is.matrix(X) && !is.list(X)) {
      stop("X parameter must be a matrix or a list.")
    }
    if (!is.matrix(Y) && !is.list(Y)) {
      stop("Y parameter must be a matrixor a list.")
    }
    if (!is.double(beta.hat) && !is.null(beta.hat)) {
      stop("The beta.hat parameter must be floating point number.")
    }
    if (!is.double(tol)) {
      stop("The tol parameter must be floating point number.")
    }
    if (is.matrix(Y)) {
      if (nrow(Y) != nrow(X)) {
        stop("Y must have the same number of rows as X.")
      }
      if (ncol(Y) < 1) {
        stop("Y must have at least 1 column.")
      }
      if (nrow(Y) != nrow(locs)) {
        stop("Y must have the same number of rows as locs.")
      }
    }
    if (is.null(scaling)) {
      if (is.matrix(locs)) {
        scaling <- rep(1, ncol(locs))
      } else {
        scaling <- rep(1, ncol(locs[[1]]))
      }
    }
    nscale <- length(unique(scaling))
    if (sum(sort(unique(scaling)) == 1:nscale) < nscale) {
      stop("scaling must consist of sequential integers between
                   1 and ncol(locs)")
    }
    if (apanasovich & nscale > 1) {
      stop("Apanasovich models require a common scale parameter")
    }
    model@scaling <- scaling
    model@nscale <- nscale
    model@apanasovich <- apanasovich
    #            if (is.matrix(locs)) {
    #                if(ncol(locs) != 2 && ncol(locs) != 3)
    # { stop("Locs must have either 2 or 3 columns.") }
    #            }
    if (is.null(covparams)) {
      model <- calc_covparams(model, locs, Y)
    }
    #            if(is.null(beta.hat)){
    #              if(is.null(ncol(Y))){
    #                ncol <- 1
    #              } else{
    #                ncol <- ncol(Y)
    #              }
    #              beta.hat <- matrix(0.0, nrow = ncol(X),ncol = ncol)
    #                if (is.matrix(X)) {
    #                    beta.hat <- matrix(0.0, nrow = ncol(X), ncol=1)
    #                }
    #                else {
    #                    beta.hat <- matrix(0.0, nrow =
    #                     ncol(superMatrix(X)), ncol=1)
    #                }
    #            }
    if (!is.double(model@covparams)) {
      stop("The covparams paramter must be a numeric vector.")
    }
    m <- model@n_neighbors
    if (m < model@min_m) {
      stop(paste("M must be at least ", model@min_m, ".", sep = ""))
    }

    if (is.list(Y)) {
      if (length(X) == 1) {
        model@X_train <- as.matrix(X[[1]])
      } else {
        model@X_train <- psych::superMatrix(X)
      }
      model@Y_train <- as.matrix(unlist(Y))
    } else {
      model@X_train <- X
      model@Y_train <- Y
    }
    if (!is.list(locs)) {
      model@locs_train <- list(locs)
    } else {
      model@locs_train <- locs
    }

    model <- specify(model, locs, m)


    if (is.null(beta.hat)) {
      beta0.glmnet <- cv.glmnet(model@X_train, model@Y_train,
        parallel = parallel,
        foldid = foldid
      )
      beta.hat <- as.matrix(predict(beta0.glmnet,
        type = "coefficients",
        s = beta0.glmnet$lambda.1se
      ))
    }
    Y.hat <- beta.hat[1, 1] + model@X_train %*% beta.hat[-1, ]
    #              dim(Y.hat) <- c(nrow(model@Y_train), ncol(model@Y_train))

    # Begining algorithm (Algorithm 1 from Messier and Katzfuss 2020)
    model@converged <- FALSE
    prev.error <- 1e10
    iter <- 1
    if (verbose) {
      cat("\n")
    }
    while (!model@converged && (iter < max_iters)) {
      model <- compute_residuals(model, model@Y_train, Y.hat)
      res_matrix <-
        matrix(model@res,
          nrow = nrow(model@Y_train),
          ncol = ncol(model@Y_train)
        )
      if (verbose) {
        cat("MSE: ", colMeans(res_matrix^2), "\n")
      }
      model <-
        estimate_theta(model, locs, optim.control, optim.method)
      # transform data to iid
      if (!model@apanasovich) {
        model <- specify(model, locs, m)
      }

      model <- transform_data(model, model@Y_train, model@X_train)
      model <- estimate_betas(model, parallel)
      min.error <- compute_error(model)
      ### Check min-error against the previous error and tolerance
      if (min.error < prev.error * tol) {
        prev.error <- min.error
        model@error <- prev.error
        beta.hat <- sparseToDenseBeta(model@linear_model)
        model@beta <- beta.hat
        # Y.hat <- as.matrix(predict(model@linear_model,
        # newx = X, s=model@linear_model$lambda[model@lambda_1se_idx]))
        Y.hat <-
          as.matrix(predict(
            model@linear_model,
            newx = model@X_train,
            s = "lambda.1se"
          ))
        covparams.iter <- model@covparams
        Vecchia.SCAD.iter <- model@linear_model
      } else {
        model@converged <- TRUE
        model@beta <- beta.hat
        model@covparams <- covparams.iter
        model@linear_model <- Vecchia.SCAD.iter
        model@error <- prev.error
      }
      if (verbose) {
        cat("\nIteration: ", iter, "\n")
        show(model)
      }
      iter <- iter + 1
    }
    return(model)
    invisible(model)
  }
)

#' estimate_betas
#'
#' Estimate the beta coefficients for a model (not called by user)
#'
#' @param model the model to estimate coeffients for
#'
#' @return A model with updated coefficients
setMethod("estimate_betas", "PrestoGPModel", function(model, parallel) {
  if (ncol(model@Y_train) > 1) {
    model@linear_model <-
      cv.glmnet(
        as.matrix(model@X_tilde),
        as.matrix(model@y_tilde),
        family = "mgaussian",
        alpha = model@alpha,
        parallel = parallel
      )
  } else {
    model@linear_model <-
      cv.glmnet(
        as.matrix(model@X_tilde),
        as.matrix(model@y_tilde),
        alpha = model@alpha,
        parallel = parallel
      )
  }
  idmin <-
    which(model@linear_model$lambda == model@linear_model$lambda.min)
  semin <-
    model@linear_model$cvm[idmin] + model@linear_model$cvsd[idmin]
  lambda_1se <-
    max(model@linear_model$lambda[model@linear_model$cvm <= semin])
  model@lambda_1se_idx <-
    which(model@linear_model$lambda == lambda_1se)
  invisible(model)
})

#' sparseToDenseBeta
#'
#' Convert the sparse beta coefficients matrix from glmnet to a dense matrix
#'
#' @param linear_model the glmnet model
#'
#' @return A dense matrix
sparseToDenseBeta <- function(linear_model) {
  coefs <- coef(linear_model)
  if (!is.list(coefs)) {
    coefs <- list(coefs)
  }
  beta_construct <-
    matrix(
      data = 0,
      nrow = coefs[[1]]@Dim[1],
      ncol = length(coefs)
    )
  # coefs[[1]]@Dim[1]+2s because dgCMatrix is 0 offset,
  # and we want to include intercept
  for (i in seq_along(coefs)) {
    for (j in seq_along(coefs[[i]]@i)) {
      k <- coefs[[i]]@i[j]
      # beta_construct[k+1,i] <- coefs[[i]]@x[j]
      beta_construct[k + 1, i] <- coefs[[i]]@x[j]
    }
  }
  # show(beta_construct)
  beta <-
    matrix(beta_construct,
      nrow = coefs[[1]]@Dim[1],
      ncol = length(coefs)
    )
  beta
}

#' compute_error
#'
#' Compute the error (log likelihood using the GP log likelihood
#' and penalty from the beta coefficients)
#'
#' @param model the PrestoGP model object
#'
#' @return The total error used in the main optimization loop
setMethod("compute_error", "PrestoGPModel", function(model) {
  ### Betas
  beta.iter <- sparseToDenseBeta(model@linear_model)
  # beta.iter <- as.numeric(coef(model@linear_model))
  # beta.iter <- do.call(cbind, coef(model@linear_model))

  ### Lambdas
  lambda.iter <- model@linear_model$lambda[model@lambda_1se_idx]

  ### Get SCAD penalty values
  # LL.vecchia.beta <- SCAD_Penalty_Loglike(beta.iter,lambda.iter)

  ### Compute log-likelihood
  # error <- model@LL_Vecchia_krig + LL.vecchia.beta[model@lambda_1se_idx[[1]]]
  error <-
    model@LL_Vecchia_krig + glmnet_penalty(beta.iter, lambda.iter, model@alpha)

  # Min error (stopping criterion) is the log-likelihood
  error
})

setMethod("transform_covariance_parameters", "PrestoGPModel", function(model) {

})
