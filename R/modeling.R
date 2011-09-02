fitNormal <- function(vals, plot=TRUE) {
  normfit <- fitdistr(vals, "normal")
  if (plot) {
    x <- seq(min(vals), max(vals), length=100)
    y <- dnorm(x, normfit$estimate[1], normfit$estimate[2])
    plot(density(vals))
    lines(x, y, col="red")
  }  
  return(normfit)
}

fitCauchy <- function(vals, plot=TRUE) {
  dfit <- fitdistr(vals, "cauchy")
  if (plot) {
    x <- seq(min(vals), max(vals), length=100)
    y <- dcauchy(x, dfit$estimate[1], dfit$estimate[2])
    plot(density(vals))
    lines(x, y, col="red")
  }
}

fitT <- function(vals, plot=TRUE) {
   dfit <- fitdistr(vals, "t")
  if (plot) {
    x <- seq(min(vals), max(vals), length=100)
    y <- dt(x, dfit$estimate[1], dfit$estimate[2])
    plot(density(vals))
    lines(x, y, col="red")
  }
}
