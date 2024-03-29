% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/datadocs.R
\docType{data}
\name{soil}
\alias{soil}
\title{Soil data from North Bavaria, Germany}
\format{
A data frame with the following columns:
\describe{
\item{x.coord:}{x coordinate (cm)}
\item{y.coord:}{y coordinate (cm)}
\item{nr:}{Sample number. The order of the sample number corresponds to the
order that the samples were collected.}
\item{moisture:}{Moisture content (Kg/Kg * 100\%)}
\item{NO3.N:}{Nitrate nitrogen (mg/Kg)}
\item{Total.N:}{Total nitrogen (mg/Kg)}
\item{NH4.N:}{Ammonium nitrogen (mg/Kg)}
\item{DOC:}{Dissolved organic carbon (mg/Kg)}
\item{N20N:}{Nitrous oxide (mg/Kg dried substance)}
}
}
\source{
The data were collected by Wolfgang Falk, Soil Physics Group,
University of Bayreuth, Germany. This data set was previously published in
the "RandomFields" R package.
}
\usage{
data(soil)
}
\description{
Soil physical and chemical data collected on a field in the Weissenstaedter
Becken, Germany
}
\details{
For technical reasons some of the data were obtained as
differences of two measurements (which are not available anymore).
Therefore, some of the data have negative values.
}
\references{
\itemize{
\item Falk, W. "Kleinskalige raeumliche Variabilitaet von Lachgas und
bodenchemischen Parameters (Small scale spatial variability of nitrous
oxide and pedo-chemical parameters)" (2000). Master thesis, University of
Bayreuth, Germany.
}
}
\author{
Martin Schlather, School of Business Informatics and Mathematics,
University of Mannheim
}
\keyword{datasets}
