plot_surface <- function(samples) {
  # library(ggplot2)
  library(scatterplot3d)
  library(rgl)

  samples <- as.data.frame(samples)
  colnames(samples) <- c("x", "y", "z")
  colors <- function(colors, values) {
    norm_values <- (values - min(values)) / diff(range(values))
    ramp <- colorRamp(colors)(norm_values)
    rgb(ramp[, 1], ramp[, 2], ramp[, 3], maxColorValue = 255)
  }
  cols <- colors(c("red", "blue"), samples$z)
  # scatterplot3d(samples[,1:3], angle=65, pch = 16, color=cols)
  # scatter3d(samples$x, samples$y, samples$z, )
  plot3d(samples$x, samples$y, samples$z, col = cols)
  # ggplot(samples, aes(x=x, y=y, z=z)) + geom_point(aes(colour=z))  + scale_colour_gradient2()
  # + scale_colour_gradientn(colours = terrain.colors(10))
  # + scale_colour_gradient2()
}
