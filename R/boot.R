library(boot)

boot.samplemean <- function(x, d, ...) {
  return(mean(x[d], na.rm=TRUE))
}

boot.median <- function(x, d) {
  return(median(x[d]))
}

boot.sum <- function(x, d) {
  return(sum(x[d]))
}

bootCI <- function(vals, stat_type="mean") {
  stat <- ""
  if (stat_type=="mean") stat <- boot.samplemean
  if (stat_type=="sum") stat <- boot.sum
  if (stat_type=="median") stat <- boot.median
  bootObj <- boot(data=vals, statistic=stat, R=100)
  bootCI <- boot.ci(bootObj, type="perc")
  percentCI <- as.array(bootCI$percent[1,4:5])
  return(percentCI)
}
