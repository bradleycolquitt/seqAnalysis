library(ggplot2)

gg_scatter <- function(data, x, y, alpha, unity_line=TRUE, density=TRUE) {
  gg <- ggplot(data, aes_string(x=x, y=y))
  gg <- gg + geom_point(alpha=I(1/alpha))
  if (unity_line) {
    gg <- gg + geom_abline(slope=1, intercept=0, linetype=2)
  }
  
  if (density) {
    breaks=seq(min(data[,x]), max(data[,x]), length.out=8)
    gg <- gg + geom_density2d(aes_string(x=x, y=y), breaks=breaks)
  }
  return(gg)
}