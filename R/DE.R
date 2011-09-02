source("~/src/R/paths.R")
library(reshape)
library(multicore)
library(foreach)


registerDoMC(cores=2)

computePvalue <- function(feature, samples) {
  data <- foreach(sample=samples, .combine="cbind") %dopar% {
    
    tmp <- read.delim(paste(feature_norm_path, 
                            paste(feature, "split", sep="_"), sample, sep="/"), header=FALSE)
    return(tmp[,4:5])
  }
  data <- data[,-3]
  colnames(data) <- c("name", "s1", "s2")
  #pvals <- daply(data, .(name), summarize, t.test(s1, s2))
  data_split <- split(data, data[,1])
  #return(data_split)
  pvals <- foreach(region=data_split, .combine="rbind") %dopar% {
  
    region_length <- nrow(region)
    if (region_length >= 20) {
    #  ind <- sample(1:region_length, 20)
    #  for (i in 2:20) {
    #    ind <- cbind(ind, sample(1:region_length, 20))
    #  }  
    #  subsamples <- lapply(c(2:3), function(column) apply(ind, 2, function(x) region[x,column]))
    #  submeans <- lapply(subsamples, function(subsample) apply(subsample, 1, mean))  
    #  p <- t.test(submeans[[1]], submeans[[2]])$p.value
    #  mean1 <- submeans[[1]]
    #  mean2 <- submeans[[2]]
    
    #region <- na.omit(region)
    ind <- region[,2] > 0 | region[,3] > 0
    #print(table(ind))
    trues <- table(ind)[2]
    if (is.na(trues)) trues <- 0
    if (trues > 10) {
      p <- t.test(region[ind,2], region[ind,3])$p.value
      mean1 <- mean(region[ind,2])
      mean2 <- mean(region[ind,3])
    } else {
      p <- 1
      mean1 <- 0
      mean2 <- 0
    }  
    return(c(as.character(region[1,1]), mean1, mean2, p))
  }
  }
  #names(pvals) <- unique(data[,1])
#  pvals <- unlist(mclapply(split(data, data[,1]), mc.cores=10, function(region) t.test(region[,1], region[,2])))
#  pvals <- tapply(data[,2:3], data[,1], function(x) ttest(x[,1], x[,2]))
  pvals <- as.data.frame(pvals)
  colnames(pvals) <- c("name", samples[1], samples[2], "ttest.p")
  pvals <- data.frame(name=pvals$name,
                     as.numeric(as.character(pvals[,samples[1]])),
                     as.numeric(as.character(pvals[,samples[2]])),
                     ttest.p=as.numeric(as.character(pvals$ttest.p)))
  colnames(pvals)[2:3] <- samples
  pvals$ttest.q <- p.adjust(pvals$ttest.p, method="BH")
  pvals$fc <- log(pvals[,2]/pvals[,3], 2)
  return(pvals)
}
