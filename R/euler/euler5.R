euler5 <- function() {
  divis <- c(20,19,18,17,16,15,14,13,12,11)
  #divis <- c(10, 9, 8, 7, 6, 5)
  #step <- prod(20,19,17,13,11)
  #step <- prod(10,7,5)
  step <- 380
  #test <- 2
  test <- 2520
  record <- 0 
  while (TRUE) {
    #a <- sapply(divis, function(x) test %% x)
    for (div in divis) {
      if (test %% div > 0) { 
        break
      } else if (div == divis[10]) {
        return(test)
      }
    }
    #if (record == 10) {
    #  return(test)
    #}
    #record <- 0 
    test <- test + step
  }
}