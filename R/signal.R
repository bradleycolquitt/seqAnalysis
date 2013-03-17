

find_min <- function(vals, span=.6, window=5) {
  # loess smooth values
  
  vals_smooth <- predict(loess(vals~c(1:length(vals)), span=span))
  # for each new value determine relationship with previous value
  # in_direction -- if new value is less than previous value, decrement in_direction (initialize to window)
  # if in_direction == 0, move to out_direction
  # if new value is greater than previous value, increment out_direction
  # if out_direction == window, mark position 'window' positions back
  
  in_direction <- window
  out_direction <- 0
  i <- 2
  while (i <= length(vals_smooth)) {
    if (in_direction == 0) {
      break
    }  
    if (vals_smooth[i] < vals_smooth[i-1]) {
      in_direction <- in_direction - 1
    } 
    i <- i + 1
  }
  
  while (i <= length(vals_smooth)) {
    if (out_direction == window) {
      break
    }  
    if (vals_smooth[i] > vals_smooth[i-1]) {
      out_direction <- out_direction + 1
    } 
    i <- i + 1
  }
  
  if (out_direction == window) {
    return(i - window)  
  } else {
    return(NA)
  }
}  