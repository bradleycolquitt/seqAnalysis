## Segments

arrow <- function(x, y, col, width) {
  lines(x, y, col=col, lwd=width)
  slope <- (y[2] - y[1]) / (x[2] - x[1])
  print(slope)
  alpha <- atan(slope)
  beta <- pi/2 - alpha
  trans.x <- cos(beta) * .5
  trans.y <- sin(beta) * .5
  xl <- x[2] - trans.x
  xt <- x[2] + cos(alpha)
  xr <- x[2] + trans.x
  yl <- y[2] + trans.y
  yr <- y[2] - trans.y
  yt <- y[2] + sin(alpha)
  print(alpha)
  print(beta)
  print(trans.x)
  print(trans.y)
  polygon(c(xl, xr), y=c(yl, yr), col=col)
}

plotSamples <- function() {
  text(12, 16, "OMP")
  text(4, 8, "Neurog1")
  text(20, 8, "ICAM1")
}

PC.setup <- function() {
  plot(1, 1, type="n", xlim=c(0,24), ylim=c(0,24))
  plotSamples()
  arrow(c(5, 10), c(9, 15), col="darkblue", width=5)
}
