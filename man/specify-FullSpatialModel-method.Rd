% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PrestoGP_Full_Spatial.R
\name{specify,FullSpatialModel-method}
\alias{specify,FullSpatialModel-method}
\title{specify}
\usage{
\S4method{specify}{FullSpatialModel}(model, locs, m)
}
\arguments{
\item{model}{A model object.}

\item{locs}{a locations matrix}

\item{m}{the number of neighbors to condition on (will always be n-1 for full models)}
}
\value{
the unmodified model
}
\description{
Do nothing, full models do not specify a conditioning set, they condition on all observations
}
