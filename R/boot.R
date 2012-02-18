library(boot)

boot.samplemean <- function(x, d, ...) {
  return(mean(x[d], na.rm=TRUE))
}

boot.median <- function(x, d) {
  return(median(x[d]))
}

bootCI <- function(vals, bound="both", stat=boot.samplemean) {
  bootObj <- boot(data=vals, statistic=stat, R=100)
  bootCI <- boot.ci(bootObj, type="perc")
  percentCI <- as.array(bootCI$percent[1,4:5])
  #if (bound=="both") {
    return(percentCI)
  #} else if (bound=="upper") {
  #  return(percentCI[2])
  #} else if (bound=="lower") {
  #  return(percentCI[1])
  #}  
}
