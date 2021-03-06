library(foreach)
library(bigmemory)
library(bigtabulate)
library(cluster) 
library(boot)
source("~/src/R/LEA/util.R")
source("~/src/R/profiles2.R")

loadROI4 <- function(fname) {
  roi.path <- "~/lib/roi"
  roi <- read.table(paste(roi.path,fname,sep="/"), sep="\t", header=FALSE,
                    colClasses=c('character', 'integer', 'integer', 'character'))
  roi[,1] <- chrVecToNum(roi[,1])
#  name_index <- c(1:nrow(roi))
#  name_code <- cbind(roi[,4], name_index)
  name_file <- paste(fname, "name_code", sep="_")
  if(!file.exists(name_file)) {
    write.table(roi[,4], file=paste(fname, "name_code", sep="_"), sep="\t", quote=FALSE,
                row.names=TRUE, col.names=FALSE)
  }
  roi[,4] <- c(1:nrow(roi))
  return(as.matrix(roi))
                   
}

mergeROIWindows <- function(roi) {
  prev.chr <- roi[1, 1]
  prev.start <- roi[1, 2]
  prev.end <- roi[1, 3]
  nrows <- length(roi[,1])
  mtx <- foreach(i = 2:nrows, .combine = rbind) %do% {
    curr.roi <- roi[i,]
    if ((curr.roi[1] == prev.chr) && (curr.roi[2] == prev.end)) {
      prev.end <- curr.roi[3]
      if (i < nrows) {
        return(NA)
      } else {
        return(c(prev.chr, prev.start, prev.end))
      }
    } else {
      to.ret <- c(prev.chr, prev.start, prev.end)
      prev.chr <- curr.roi[1]
      prev.start <- curr.roi[2]
      prev.end <- curr.roi[3]
      return(to.ret)
    }
  }
  colnames(mtx) <- c("chr", "start", "end")
  rownames(mtx) <- 1:length(mtx[,1])
  return(mtx[which(!is.na(mtx[,1])),])
}

diffEnrichmentROI <- function(enrich1, enrich2) {
  roi1.sum <- sum(enrich1)
  roi2.sum <- sum(enrich2)
  p.val <- foreach(i=icount(length(enrich1)), .combine=c) %do% {    
    cont.table <- matrix(c(enrich1[i], roi1.sum - enrich1[i],
                           enrich2[i], roi2.sum - enrich2[i]), nrow=2, ncol=2)
    return(fast.fisher(cont.table)$p.val)
  }
  return(p.val)
}

prepareMPforCluster <- function(roi=NULL, data_type="ams_A", N=2000, write=FALSE) {
  vals <- extractMPvalues(roi, data_type, write)
  #return(vals)
  vals <- removeNA(vals)
  vals.filt <- filterClusteredGenes(vals)  
  vals.filt <- vals.filt[apply(vals.filt, 1, prod) > 0, ]                                        
  if (N > 0) {vals.filt <- threshMPbyVar(vals.filt, N)}
  return(vals.filt)
}

filterClusteredGenes <- function(vals) {
  clust <- scan(file="~/lib/roi/clusteredGenes_root", what=character())
  index <- unlist(sapply(clust, function(x) grep(x, rownames(vals))))
  vals.filt <- vals[-index,]
  vals.filt <- vals.filt[-grep("Rik$|Rik[0-9]", rownames(vals.filt)),]
  return(vals.filt)
}

collectClusteredGenes <- function(roi) {
  genes <- read.delim("~/lib/roi/refgene", header=FALSE, as.is=T)
  names <- genes[,4]
  names.split <- sapply(names, strsplit, "[0-9]{1,4}$|ps1$")
  names.root <- unlist(lapply(names.split, function(x) x[1]))
  names.table <- table(names.root)
  names.clust <- names(names.table[names.table>=10])
  return(names.clust)
  names.clust.ex <- names(unlist(sapply(names.clust, function(x) grep(x, names))))
  return(names.clust.ex)
                                        #  return(names%in%names.clust.ex)
  names.rik <- names[grep("Rik$|Rik[0-9]", names)]
  genes <- genes[!names%in%names.rik,]
  #return(names.clust.ex)
  #return(genes[,4]%in%names.clust.ex)
  #return(genes[,4])
  genes.noclust <- genes[!genes[,4]%in%names.clust.ex,]
  return(genes.noclust)
}

clusterByClara <- function(vals, k_range=c(2,10), plot=TRUE) {
  claras <- list()
  range <- seq(k_range[1], k_range[2])
  for(k in range) {
    cl <- clara(vals, k=k)
    claras <- c(claras, list(cl))
  }
  if (plot) {
    x11()
    x <- ceiling(length(claras)/3)
    print(x)
    par(mfrow=c(x, 3))
    for(index in 1:length(claras)) {
      clusplot(claras[[index]])
    }
  }
 
 return(claras)
}

clusterByHclust <- function(vals, dist="euclidean", method="ward") {
  vals_dist <- dist(vals, dist)
  h <- hclust(vals_dist, method=method)
  return(h)
}

clusterByHopach <- function(vals=NULL) {
  vals.dist <- distancematrix(vals, "cosangle")
  vals.hobj <- hopach(vals, dmat = vals.dist)
  return(list(vals.dist, vals.hobj))
}

clusterByMona <- function(bt=NULL, key=NULL, fname=NULL) {
  bt.mona <- mona(bt)
  step <- bt.mona$step
  order.lab <- bt.mona$order.lab
  prev <- 0
  clusters <- list()
  for(i in 1:length(step)) {
    if (step[i] == key) {
      clusters <- c(clusters, list(order.lab[(prev+1):i]))
      prev <- i
    }
  }
  if (!is.null(fname)) {
    file_path <- "~/analysis/binary/mona"
    for(i in 1:length(clusters)) {
      write.table(clusters[[i]], file=paste(file_path, paste(fname, "step", key, "cluster", i, sep="_" ), sep="/"),
                  quote=FALSE, row.names=FALSE, col.names=FALSE)
    }
  }  
  return(clusters)
}

prComp <- function(vals) {
  result <- prcomp(vals, retx=TRUE, center=TRUE, scale.=TRUE)
  return(result)
}
extractMPvalues <- function(roi=NULL, data_type="ams_A", write=FALSE) {
  data_path <- "~/analysis/mprofiles/"
  combined_names <- c(hmedips, medips)
  combined <- lapply(combined_names, function(x) {
    d <- read.delim(paste(data_path, paste(x, roi, sep="_"), sep=""), header=TRUE)
    return(d[,data_type])
  })
  roi_file <- loadROI4(roi)
  roi_file[,4] <- decodeNames("refgene_roi4_name_code", roi_file[,4])
  combined <- matrix(unlist(combined), c(length(combined[[1]]), length(combined_names)),
                     byrow=FALSE, dimnames=list(roi_file[,4], combined_names))
  return(combined)
}

threshMPbyVar <- function(vals, N) {
  vars <- apply(vals, 1, var)
  q <- quantile(vars, (nrow(vals) - N) / nrow(vals))
  #return(vars[vars >= q])
  vals.selected <- vals[vars >= q, ]
  return(vals.selected)
}

removeNA <- function(mat) {
  mat <- mat[apply(mat, 1, NAtest),]
  return(mat)
}

NAtest <- function(x) {
  bool <- !is.na(x) && !is.nan(x) && is.finite(x)
  return(bool)
}


threshMedipsByVal <- function(roi=NULL, data_type="ams_A", thresh=700, write=FALSE, binary=TRUE) {
  data_path <- "~/analysis/mprofiles/"
  combined_names<- c(hmedips, medips)
  combined <- lapply(combined_names, function(x) {
    d <- read.delim(paste(data_path, paste(x, roi, sep="_"), sep=""), header=TRUE)
    if (thresh > 0) {
      d <- d[is.finite(d[,data_type]),]
      d <- d[d[,data_type] >= thresh, ]}
    return(d)
  })
  names(combined) <- combined_names
  return(combined)
  roi_file <- loadROI4(roi)
  roi_file[,4] <- decodeNames("refgene_roi4_name_code", roi_file[,4])

  if (thresh == 0 && !binary) {
    combined <- foreach (val=combined, .combine="cbind") %do% {return(val[,data_type])}
    rownames(combined) <- roi_file[,4]
    colnames(combined) <- combined_names
    if (write) {
      out_path <- paste("~/analysis/mprofiles/roi/", paste(roi, data_type, sep="_"), sep="")
      write.table(output, file=out_path, quote=FALSE, sep="\t")
    }
    return(combined)
  }
  output <- matrix(0, nrow=nrow(roi_file), ncol=length(combined_names), dimnames=list(roi_file[,4]))
  for (index in 1:length(combined)) {
    matches <- match(rownames(output), rownames(combined[[index]]))
    output[matches, index] <- 1
  }
  colnames(output) <- combined_names
  if (write) {
    out_path <- paste("~/analysis/binary/", paste(roi, data_type, thresh, sep="_"), sep="")
    write.table(output, file=out_path, quote=FALSE, sep="\t")
  }
  return(output)
}

threshRNAbyVal <- function(thresh=0) {
  rna_path <- "~/data/rna"
  rna_files <- c("ompgfp.txt", "ngnhigh.txt")
  filtered_rna <- lapply(rna_files, function(file) {
    d <- read.delim(paste(rna_path, file, sep="/"), header=FALSE, as.is=TRUE)
    d[,2] <- log(d[,2], 2)
    d <- d[d[,2] >= thresh,]
    return(d)
  })
  names(filtered_rna) <- c("omp", "ngn")
  return(filtered_rna)
}

makeBinaryTable <- function(vals, reference, code=TRUE) {
  ref <- reference
  if (code) {
    ref <- loadROI4(reference)
    ref[,4] <- decodeNames(paste(reference, "roi4_name_code",sep="_"), ref[,4])
                         }
  output <- matrix(0, nrow=nrow(ref), ncol=length(vals),
                   dimnames=list(ref[,4], names(vals)))
  for (index in 1:length(vals)) {
    matches <- match(rownames(output), vals[[index]])
    output[matches, index] <- 1
  }
  zero_index <- apply(output, 1, sum)
  output <- output[zero_index>0,]
  output <- filterClusteredGenes(output)
  return(output)
}


makeBinaryTable.HS <- function(vals, reference, code=TRUE) {
  ref <- reference
  if (code) {
    ref <- loadROI4(reference)
    ref[,4] <- decodeNames(paste(reference, "roi4_name_code",sep="_"), ref[,4])
                         }
  output <- matrix(0, nrow=nrow(ref), ncol=length(vals),
                   dimnames=list(ref[,4], names(vals)))
  for (index in 1:length(vals)) {
    matches <- match(rownames(output), vals[[index]][,4])
    output[vals[[index]][,4], index] <- 1
  }
  zero_index <- apply(output, 1, sum)
  output <- output[zero_index>0,]
  #output <- filterClusteredGenes(output)
  return(output)
}

binaryHS <- function(path) {
  dir_files <- list.files(path)
  sample_files <- dir_files[grep("inter", dir_files)]
  ref <- dir_files[grep("DNAse.txt", dir_files)]
  file_vals <- lapply(sample_files, function(x) {
    d <- read.delim(paste(path,x,sep="/"), header=F)
    return(d)
  })
  sample_names <- sapply(sample_files, strsplit, "DNAse_inter_")
  sample_names <- unlist(lapply(sample_names, function(x) x[[2]]))
  names(file_vals) <- sample_names
  ref_file <- read.delim(paste(path,ref,sep="/"), header=F)
  bt <- makeBinaryTable.HS(file_vals, ref_file, code=FALSE)
}

binaryDNAmodAndRNA <- function(roi=NULL, thresh.medips=700, thresh.rna=0) {
  medips <- threshMedipsByVal(roi=roi, thresh=thresh.medips)
  medips.names <- lapply(medips, rownames)
  rna <- threshRNAbyVal(thresh=thresh.rna)
  rna.names <- lapply(rna, function(x) x[,1])
  bt <- makeBinaryTable(c(medips.names, rna.names), roi)
  return(bt)
}
 
calcBTdist <- function(bt, distFunc=XOR, asMatrix=FALSE) {
  comb <- combn(colnames(bt), 2)

  XOR <- function(x, y) {
    z <- xor(x,y)
    score <- (length(x) - sum(z)) / length(x)
    return(score)
  }

  distance <- function(x, y) {
    val <- dist(rbind(x, y))
    return(val)
  }
  
  scores <- apply(comb, 2, function(x) {
    val <- do.call(XOR, list(bt[,x[1]], bt[,x[2]]))
    return(val)
  })
  
  out_names <- apply(comb, 2, paste, collapse=" to ")
  names(scores) <- out_names
  if (asMatrix) {
    mat <- fillMatrix(colnames(bt),scores)
    out_names <- colnames(bt)
    dimnames(mat) <- list(out_names, out_names)
    return(mat)
  }
  return(scores)
}

fillMatrix <- function(set, vals) {
   mat <- matrix(1, nrow=length(set), ncol=length(set))
   .filler <- function(mat) {
     prev <- 0
     for(i in 1:(length(set) - 1)) {
       step <- prev + length(set) - i
       mat[i, (i+1):ncol(mat)] <- vals[(prev + 1):step]
       prev <- step
    }
    return(mat)
   }
   mat <- .filler(mat)
   mat <- .filler(t(mat))
   return(mat)
   
  }



testKmeans <- function(mydata) {
  wss <- (nrow(mydata)-1)*sum(apply(mydata,2,var))
  for (i in 2:15) wss[i] <- sum(kmeans(mydata,
                                      centers=i)$withinss)
  plot(1:15, wss, type="b", xlab="Number of Clusters",
       ylab="Within groups sum of squares")
}

##Not currently used
binaryTable <- function(files, roi.name) {
  data.path <- "~/analysis/binary"
  files <- list.files(data.path)
  roi <- loadROI4(roi.name)
  output <- matrix(0, nrow = length(roi), ncol = length(files))
  for (index in 1:length(files)) {
    scores <- read.delim(paste(data.path, files[index], sep="/"), header=FALSE,
                         colClasses = c('numeric', 'numeric'))
    output[scores[,1], index] <- 1
  }
  colnames(output) <- files
  return(output)     
  
}

clusterImageByPos <- function(data, pos, centers) {
  data.k <- kmeans(data[,pos], centers=centers)
  data.means <- lapply(split(data, data.k$cluster), function(x) apply(x, 2, mean, na.rm=T))
  return(data.means)
}

clusterImage.samplePivot <- function(set, anno, pos, centers, pivot) {
  if (set == "cells") {
    samples <- c(hmedips, medips)
  }
  data <- lapply(samples, function(x) read.delim(paste(profile2.path, x, "images", anno, sep="/"), header=FALSE, row.names=1))
  names(data) <- samples
  pivot.k <- kmeans(data[[pivot]][,pos], centers=centers)
  means <- lapply(data, function(datum) {
    sample.mean <- lapply(split(datum, pivot.k$cluster), function(x) apply(x, 2, mean, na.rm=T))
  })
  return(means)
}

clusterImage.plot <- function(samples.means) {
  for(i in 1:length(samples.means)) {
    plot(1,1, type="n", xlim=c(1,length(samples.means[[i]][[1]])), ylim=c(0,1))
    curr.sample <- samples.means[[i]]
    for (j in 1:length(curr.sample)) {
      lines(curr.sample[[j]], col=j)
    }
  }
}

## Measure of similarity between hierarchical clusterings
## From Fowlkes and Mallows, JASA, 1983

testClusterSim <- function(cluster1, cluster2, range=c(2,10)) {
  B_ks <- foreach (k=range) %do% {
  k <- 2
    #cluster1 <- cutree(hc1, k)
    #cluster2 <- cutree(hc2, k)
    
    simmatrix <- matrix(nrow=k, ncol=k)
    for (i in 1:k) {
      names1 <- names(cluster1[cluster1==i])
      for (j in 1:k) {
        names2 <- names(cluster2[cluster2==j])
        simmatrix[i,j] <- length(intersect(names1, names2))
      }
    }

    m_i <- apply(simmatrix, 1, sum)
    m_j <- apply(simmatrix, 2, sum)
    n <- sum(simmatrix)
    P_k <- sum(m_i ^ 2) - n
    Q_k <- sum(m_j ^ 2) - n
    T_k <- sum(simmatrix ^ 2) - n
    B_k <- T_k / sqrt(P_k * Q_k)
    return(B_k)
 }
  #return(B_k)
  names(B_ks) <- range
  return(B_ks)

}

## Cut trees into k clusters
## If permuting, test range of permutation fractions for a given k number of clusters
testClusterSim.range <- function(hc1, hc2, range=seq(2, 10), perm=c(.05, .1, .2, .5, .8, 1)) {
  B_ks <- foreach(k=range) %do% {
    cluster1 <- cutree(hc1, k)
    cluster2 <- cutree(hc2, k)
    out <- c()
    if (length(perm) > 1) {
      for (fraction in perm) {
        cluster2 <- sampleSubset(cluster1, fraction)
        out <- c(out, testClusterSim(cluster1, cluster2, k))
      }  
    } else {
      out <- testClusterSim(cluster1, cluster2, k)
    }
    return(out)
  }
  names(B_ks) <- range
  return(B_ks)
}


testClusterSim.bootInterface <- function(hc1, hc2, k=2, R=10) {
  cluster1 <- cutree(hc1, k)
  cluster2 <- cutree(hc2, k)
  #hc_boot <- boot(data=cluster2, statistic=c, R=R)
  #hc_boot.array < boot.array(hc_boot)
  #out <- foreach (row=hc_boot.array) %do% {
  i <- 1
  out <- vector(length=R)
  while (i <= R) {
    cluster2 <- cluster2[sample(1:length(cluster2), length(cluster2))]
    out[i] <- testClusterSim(cluster1, cluster2, k=k)
    
    i <- i + 1
  }  
  #}
  return(out)
}

sampleSubset <- function(data, fraction) {
  num_sample <- round(length(data) * fraction)
  #print(num_sample)
  sample_ind <- sample(1:length(data), num_sample)
  #print(sample_ind)
  sample_ind_new <- as.numeric(sample(as.character(sample_ind), length(sample_ind)))
  #print(sample_ind_new)
  data[sample_ind_new] <- data[sample_ind]
  return(data)
}

computeClusterStats.hc <- function(d, hc, range=c(2:10)) {
  stat_list <- list()
  for (i in range) {
    print(i)
    ct <- cutree(hc, i)
    stat_list <- c(stat_list, list(cluster.stats(d, ct)))
  }
  return(stat_list)
}

hmedips <- paste(c("omp","ngn","icam"),"hmedip.bed",sep="_")
medips <- paste(c("omp","ngn","icam"),"medip.bed",sep="_")
rna <- paste(c("omp","ngn","icam"),"rna",sep="_")
