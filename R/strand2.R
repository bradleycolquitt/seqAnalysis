library(plyr)

# Compute mean values for flanking and gene bodies 
compute_mean_gene <- function(data) {
  mid <- apply(data[,51:100], 1, mean, na.rm=TRUE)
  flank <- apply(data[,c(1:50,101:150)], 1, mean, na.rm=TRUE)
  df <- data.frame(pos=rep(c("Body","Flank"), each=length(mid)), val=c(mid, flank))
  return(df)
}

# Compute mean values for flanking and gene bodies for set of matrices
compute_mean_gene_set <- function(data_list) {
  dfs <- lapply(data_list, compute_mean_gene)
  dfs_un <- ldply(dfs)
  dfs_un[,1] <- as.factor(as.numeric(dfs_un[,1]))
  colnames(dfs_un)[1] <- "group"
  #nrows <- unlist(lapply(dfs, function(x) return(dim(x)[1])))
  #print(names(data_list))
  #dfs_un$group <- rep(names(data_list), times=nrows)
  return(dfs_un)
}