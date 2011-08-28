
## Read in intersections and make count matrix

intersectToMatrix <- function() {
  file.path <- "~/data/ecmp/reads_count_inter_MACS_DNAse"
  files <- list.files(file.path)
  files <- files[grep("norm_wg", files)]

  path <- paste(file.path, files[1], sep="/")
  counts <- readInIntersections(path) 
  for(file in files[2:length(files)]) {
    path <- paste(file.path, file, sep="/")
    counts <- cbind(counts, readInIntersections(path))
  }
  file.names <- sapply(files, function(x) strsplit(x, "_reads"))
  file.names <- unlist(lapply(file.names, function(x) x[1]))
  colnames(counts) <- file.names
  d <- read.delim(path, header=FALSE)
  rownames(counts) <- d[,4]
  return(counts)
}

readInIntersections <- function(file) {
  data <- read.delim(file, header=FALSE)
  return(data[,6])
}

clusterByPAM <- function(data) {
  pam.results <- lapply(c(2:7), function(x) pam(data, x))
  png(file="/media/XTREME!/intersections/hs_pam_clusters.png")
  par(mfrow=c(2,3))
  for(i in 1:length(pam.results)) plot(pam.results[[i]], which.plots=1)
  dev.off()
  clusters <- lapply(pam.results, function(x) x$clustering)
  for(i in 1:length(clusters)) {
    write.table(clusters[[i]], file=paste("/media/XTREME!/intersections/hs_pam_cluster_", i +1, sep=""), quote=FALSE, sep="\t", row.names=TRUE, col.names=FALSE)
  }
}

clusterByHclust <- function(data) {
  data.dist <- dist(data)
  data.hclust <- hclust(data.dist)
  png(file="/media/XTREME!/intersections/hs_hclust.png")
  par(mfrow=c(1,1))
  plot(data.hclust)
  dev.off()
  return(data.hclust)
   }
