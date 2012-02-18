library(MASS)

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

getNormQ <- function(fit.param, q) {
  d <- qnorm(q, mean=fit.param$estimate[1], fit.param$estimate[2])
  return(d)
}

permutationTest <- function(x, y, N=1000) {
  x_y <- c(x, y)
  t_obs <- abs(mean(x, na.rm=TRUE) - mean(y, na.rm=TRUE))
  t_perm <- vector("numeric", N)
  
  for (i in 1:N) {
    ind <- sample(1:length(x_y), length(x_y))
    #x_y_perm <- x_y[ind]
    half_1 <- 1:(length(x_y)/2)
    #print(half_1)
    half_2 <- (length(x_y)/2 + 1):length(x_y)
    #print(half_2)
    t_perm[i] <- abs(mean(x_y[ind[half_1]], na.rm=TRUE) - mean(x_y[ind[half_2]], na.rm=TRUE))
    #stop
  }
  #print(t_obs)
  #return(t_perm)
  above <- table(t_perm >= t_obs)
  if (length(above) == 2) {
    p_value <- above[2]/N
  } else {
    p_value <- 0
  }
  return(p_value)
  print(table(t_perm >= t_obs))
  #return(t_perm)
  #print(above)
  #p_value <- above/N
  #return(p_value)
}

## Accepts two matrices, tests pairwise by column
permutationTest.mat <- function(d1, d2, chunkSize=5, N=10000, adjust="BH") {
  registerDoMC(cores=10)
  ncols <- ncol(d1)
  i <- 0
  out <- foreach(column=isplitIndices(ncols, chunkSize=chunkSize), .combine="c") %dopar% {
    if (chunkSize > 1) {
      x <- apply(d1[,column], 1, mean, na.rm=TRUE)
      y <- apply(d2[,column], 1, mean, na.rm=TRUE)
    } else {
      x <- d1[,column]
      y <- d2[,column]
    }
    permutationTest(x, y, N=N)
  }
  if (!is.null(adjust)) p.adjust(out, method="BH")
  return(out)
}

## Accepts list of matrices and performs a pairwise comparison on all combinations
## Returns list of column pvalues for each comparison
permutationTest.matList <- function(data, chunkSize=5, N=10000, adjust="BH") {
  comparisons <- combn(names(data), 2)
  comparisons_labels <- apply(comparisons, 2, paste, collapse="_")
  out <- apply(comparisons, 2, function(comparison) {
    print(comparison)
    permutationTest.mat(data[[comparison[1]]], data[[comparison[2]]], chunkSize=chunkSize, N=N)
  })
  colnames(out) <- comparisons_labels
  return(out)
}

plotSigLine <- function(sig_values, thresh=0.05, step, yval, ...) {
  sig_positions <- which(sig_values <= thresh)
  if (length(sig_positions > 0)) {
    for(i in 1:length(sig_positions)) {
      lines(c(sig_positions[i] * step, sig_positions[i] * step + step ), rep(yval, times=2), lwd=2, ...)
    }
  } else {
    print("No significant positions.")
  }  
}
