source("~/src/R/paths.R")
library(reshape)
library(multicore)
library(foreach)


registerDoMC(cores=10)

computePvalue <- function(set, feature, samples) {
  data <- foreach(sample=samples, .combine="cbind") %dopar% {
    
    tmp <- read.delim(paste(feature_norm_path, 
                            paste(feature, "split", sep="_"), sample, sep="/"), header=FALSE)
    return(tmp[,4:5])
  }
  data <- data[,-3]
  colnames(data) <- c("name", "s1", "s2")
  #pvals <- daply(data, .(name), summarize, t.test(s1, s2))
  data_split <- split(data, data[,1])
  pvals <- foreach(region=data_split, .combine="cbind") %dopar% {
    p <- t.test(region[,2], region[,3])$p.value
    mean1 <- mean(region[,2])
    mean2 <- mean(region[,3])
    return(c(region[1,1], mean1, mean2, p))
  }
  #names(pvals) <- unique(data[,1])
#  pvals <- unlist(mclapply(split(data, data[,1]), mc.cores=10, function(region) t.test(region[,1], region[,2])))
#  pvals <- tapply(data[,2:3], data[,1], function(x) ttest(x[,1], x[,2]))
  pvals[,4] <- p.adjust(pvals[,4], method="BH")
  return(pvals)
}
