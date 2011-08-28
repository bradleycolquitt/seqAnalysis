library(plyr)
library(ggplot2)
library(splines)
deseas <- function(var, month) {
  val <- resid(lm(var~factor(month))) + mean(var, na.rm=TRUE)
  return(val)
}
