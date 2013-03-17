# install.packages("BB") # if you already don't have it.
require(BB)

# generate test data to check model efficiency
set.seed(100) 
gamma.nb.mix <- function(n, prob=0.5, size1=10, mu1=0, size2=20, mu2=4) {
  u <- runif(n) 
  out <- apply( as.matrix(u), 1, function(x) ifelse(x<=prob, rnbinom(1, size=size1, mu=mu1), rnbinom(1, size=size2, mu=mu2) ) )
  out 
}
out <- gamma.nb.mix(1000)

# maximizing function
fit_mix_nbinom <- function(data) {
f = function(size1, mu1, size2, mu2, prob) { 
  -sum(log(prob*dnbinom(data, size=size1, mu=mu1) + 
    (1-prob)*dnbinom(data, size=size2, mu=mu2))) 
}

# arbitrary start (probability parameter is important to remain as close as possible)
start0 <- list("size1"=20, "mu1"=4, "size2"=20, "mu2"=4, "prob"=0.5) 

fcn <- function(x) f(x[1], x[2], x[3], x[4], x[5])   
res <- optim(par=as.numeric(start0), fn=fcn, control=list( maxit=10000 ) , method="BFGS" )
return(res)
}