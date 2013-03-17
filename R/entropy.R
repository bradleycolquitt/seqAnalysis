discretize <- function(vals, breaks="FD") {
  if (is.character(breaks)) {
    breaks <- hist(vals, breaks=breaks, plot=FALSE)$breaks
  } 
  cut(vals, breaks=breaks, labels=FALSE, include.lowest=TRUE)
  
}