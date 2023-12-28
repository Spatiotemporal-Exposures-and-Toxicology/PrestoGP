library(fields)

quilt_plot <- function(locs, values) {
  n <- length(values)
  fields::quilt.plot(locs, values, nx = sqrt(n), ny = sqrt(n))
}
