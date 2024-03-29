% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PrestoGP_Multivariate_Vecchia.R
\docType{class}
\name{MultivariateVecchiaModel-class}
\alias{MultivariateVecchiaModel-class}
\alias{MultivariateVecchiaModel}
\title{Multivariate Vecchia PrestoGP model class}
\description{
This class is used to create multivariate models with a likelihood function
conditioned on a subset of the observations (i.e., Vecchia models). See
\code{\link{PrestoGPModel}} for a description of the slots.
}
\examples{
pgp.mmodel <- new("MultivariateVecchiaModel", n_neighbors = 25)
}
\seealso{
\code{\link{PrestoGPModel}}
}
