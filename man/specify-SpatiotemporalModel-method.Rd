% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PrestoGP_Vecchia_Spatiotemporal.R
\name{specify,SpatiotemporalModel-method}
\alias{specify,SpatiotemporalModel-method}
\title{specify}
\usage{
\S4method{specify}{SpatiotemporalModel}(model, locs, m)
}
\arguments{
\item{model}{The model to specify}

\item{locs}{the locations matrix}

\item{m}{the number of neighbors to condition on}
}
\value{
a model with a specified conditioning set
}
\description{
Specify the conditioning set using m nearest neighbors.
}
