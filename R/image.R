library(itertools)
library(gplots)
library(plyr)
library(foreach)
library(doMC)
library(RColorBrewer)

source("~/src/seqAnalysis/R/paths.R")
source("~/src/seqAnalysis/R/plotUtil.R")

ybb <- colorRampPalette(c("yellow", "black", "blue"), space="rgb")

## Generates matrix of samples values with annotation observations as rows and
##    positions as columns
##    Arguments:  data - profile data with [chr start stop name position strand value]
##                data_type - directory and normalization type
positionMatrixRep <- function(data, group2, data_type="raw") {
  matrices <- lapply(data, function(d) {  
    groups <- unique(d[,5])
    names <- unique(d[,4])
    print("Splitting...")
    data_groups = split(d, d[,5])
    print("Combining...")
    out <- foreach(ind=icount(length(data_groups)), .combine="cbind") %dopar% {
      data_group <- data_groups[[ind]]
      m <- match(names, data_group[,4])
      return(data_group[m, 7])
    }
    rownames(out) <- names
    return(out[!is.na(rownames(out)),])
  })
  #return(matrices)
  if (length(matrices) > 1) {
    return(matrix_el_sum(matrices) / length(matrices))
  } else {
    return(matrices[[1]])
    
  }
}

## Generates matrix of samples values with annotation observations as rows and
##    positions as columns
##    Arguments:  data - profile data with [chr start stop name position strand value]
##                data_type - directory and normalization type
positionMatrix <- function(data, group2, data_type="raw") {
  groups <- unique(data[,5])
  names <- unique(data[,4])
  print("Splitting...")
  data_groups = split(data, data[,5])
  print("Combining...")
  out <- foreach(ind=icount(length(data_groups)), .combine="cbind") %dopar% {
    data_group <- data_groups[[ind]]
    m <- match(names, data_group[,4])
    return(data_group[m, 7])
  }
  rownames(out) <- names
  return(out[!is.na(rownames(out)),])
}

matrix_el_sum <- function(x) {
  out <- x[[1]]
  for (i in 1:length(x)) {
    out <- out + x[[i]]
  }
  return(out)
}

positionMatrix.group <- function(anno, sample, data_type="rpkm/mean", group2=NULL, rm.outliers=0) {
  data_path <- paste(profile2.path, "norm", data_type, anno, sep="/")
  #data <- read.delim(paste(data_path, sample, sep="/"), header=FALSE, colClasses=image.classes)
  
  group2_data <- read.delim(paste(group2.path, group2, sep="/"), header=FALSE)
  
  out_path <- paste(profile2.path, "norm", data_type, anno, "images", sep="/")
  out_name <- paste(sample, group2, sep="_")
  if (!file.exists(out_path)) dir.create(out_path)
  
  if (!file.exists(paste(out_path, paste(out_name, group2_data[1,2], sep="_"), sep="/"))) {
    print(sample)
    data_path <- paste(profile2.path, "norm", data_type, anno, sep="/")
    data <- read.delim(paste(data_path, sample, sep="/"), header=FALSE, colClasses=image.classes)
    #return(data)
    group2_data <- read.delim(paste(group2.path, group2, sep="/"), header=FALSE)
    
    if (rm.outliers > 0) {
      thresh <- quantile(data[,7], probs=c(rm.outliers, 1-rm.outliers), na.rm=TRUE)
      data <- data[data[,7] >= thresh[1] & data[,7] <= thresh[2], ]
    }
    #return(match(as.character(data[,4]), as.character(group2_data[,1])))
    group2_ex <- group2_data[match(data[,4], as.character(group2_data[,1])),]
    #return(group2_ex)
    data_split <- split(data, group2_ex[,2])
    #return(data_split)
    a <- laply(names(data_split), function(x) {
      pos_matrix <- positionMatrix(data_split[[x]], data_type=data_type)
      write.table(pos_matrix, file=paste(out_path, paste(out_name, x, sep="_"), sep="/"),
                  quote=FALSE, sep="\t", col.names=FALSE)
      a <- gc()
    }, .progress="text")
  }  
}
  
## Run positionMatrix for specified sample set and given annotation
##     Arguments:  anno - annotation file
##                 data_type - directory and normalization type
positionMatrixRep.all <- function(anno, data_type="unnorm/mean", rep=FALSE) {
  
  data_path <- paste(profile2.path, "norm", data_type, anno, sep="/")
  samples <- list.files(data_path)
  
  if (rep) {
    samples <- group_by_rep(samples)
    #print(samples)
  }
  #return(samples)
  for (i in 1:length(samples)) {
    
    sample_name <- names(samples)[i]
    sample <- samples[[i]][[1]]
    print(sample_name)
    print(sample)
    if (!(sample_name == "images" | sample_name == "profiles")) {
      out_path <- paste(profile2.path, "norm", data_type, anno, "images", sep="/")
      if (!file.exists(paste(out_path, sample_name, sep="/"))) {
        if (length(sample) == 1) {
          data <- read.delim(paste(profile2.path, "norm", data_type, anno, sample, sep="/"),
                       header=FALSE, colClasses=image.classes)
        } else {
          data <- lapply(sample, function(s) read.delim(paste(profile2.path, "norm", data_type, anno, s, sep="/"),
                                                        header=FALSE, colClasses=image.classes))  
        }  
        pos_matrix <- positionMatrix(data, data_type=data_type)
 
        if (!file.exists(out_path)) dir.create(out_path)
        write.table(pos_matrix, file=paste(out_path, sample_name, sep="/"),
                quote=FALSE, sep="\t", col.names=FALSE)
        a <- gc()
      }
    }
  }  
}

## Run positionMatrix for specified sample set and given annotation
##     Arguments:  anno - annotation file
##                 data_type - directory and normalization type
positionMatrix.all <- function(anno, data_type="unnorm/mean") {
  
  data_path <- paste(profile2.path, "norm", data_type, anno, sep="/")
  samples <- list.files(data_path)
  
  for (sample in samples) {
    if (sample != "images" & sample != "profiles") {
      out_path <- paste(profile2.path, "norm", data_type, anno, "images", sep="/")
      if (!file.exists(paste(out_path, sample, sep="/"))) {
        print(sample)
        data <- read.delim(paste(profile2.path, "norm", data_type, anno, sample, sep="/"),
                         header=FALSE, colClasses=image.classes)
        pos_matrix <- positionMatrix(data, data_type=data_type)
          
        if (!file.exists(out_path)) dir.create(out_path)
        
        write.table(pos_matrix, file=paste(out_path, sample, sep="/"),
                  quote=FALSE, sep="\t", col.names=FALSE)
        
        a <- gc()
      }
    }
  }  
}


makeImage <- function(sample, anno, data_type="unnorm", image=FALSE, heat=TRUE, range=NULL, average=NULL) {
  image.path <- paste(profile2.path, "norm", data_type, anno, "images", sample, sep="/")
  print(image.path)
  data <- as.matrix(read.delim(image.path, header=FALSE, row.names=1, as.is=T))
  #mp.data <- prepMP(mp.data)
  #mp.data <- thresh
  #out <- MP.positionMatrix(mp.data)
  if (image) {
    if (heat) {
      MP.heat(data, range, average)
    } else {
      MP.image(data)
    }  
  }
  n <- which(is.na(data), arr.ind=T)
  data[n] <- 0
  return(data)
}


makeImage.all <- function(set="d3a", anno, data_type="unnorm/mean") {
  if (set=="d3a") {
    samples <- samples.d3a
  } else if (set=="cells") {
    samples <- samples.cells_norm
  } else if (set=="cells_rpkm") {
    samples <- list("omp_hmc_120424_rpkm", "ngn_hmc_rpkm", "icam_hmc_rpkm",
                    "omp_mc_rpkm", "ngn_mc_rpkm", "icam_mc_rpkm")
  }
  
  out <- lapply(samples, function(sample) {
    makeImage(sample, anno, data_type=data_type, image=FALSE)
  }) 
  names(out) <- samples
  n <- which(is.na(out), arr.ind)
  out[n] <- 0
  return(out)
  
}

MP.processImage <- function(vals) {
  v <- apply(vals, 2, pseudoCountNorm)
  #v <- na.omit(v)
  v <- sqrt(v)
  return(v)
}
MP.image <- function(vals) {
  X11()
  image(t(vals),
        #col=rev(heat_hcl(50, c=c(80,30), l=c(30,90), power=c(1/5, 1.3))))
        col=greenred(50))
  
}

MP.heat <- function(data, density="none", range=NULL, average=NULL, fname=NULL, color="blues") {
  require(gplots)
  if (is.null(fname)) {
    #x11()
  } else {
    pdf(file=paste("~/s2/analysis/profiles/plots", fname, sep="/"), 12,12)
    #pdf(file=fname, 12, 12)
  }
  
  if (is.null(range)) {
    print(min(data, na.rm=T))
    print(max(data, na.rm=T))
    breaks <- seq(min(data, na.rm=T), max(data, na.rm=T), length.out=101)
  } else {
    if (length(range)==2) {
      breaks <- seq(range[1], range[2], length.out=101)
    } else if (length(range==4)) {
      breaks <- c(range[1], seq(range[2], range[3], length.out=99), range[4])
    }  
  }

  if (!is.null(average)) {
    #ind <- rep(seq(1, nrow(data)/average), each=average)
    data <- foreach(group=isplitRows(data, chunkSize=average), .combine="rbind", .inorder=TRUE) %dopar% {
      apply(group, 2, mean, trim=.2, na.rm=TRUE)
    }
  }
  #return(data)
  lab.palette <- c()
  if (color=="blues") {
    lab.palette <- colorRampPalette(brewer.pal(5,"Blues"), space="Lab")
  } else if (color=="ybb") {
    lab.palette <- ybb
  }  
  heatmap.2(data, Rowv=FALSE, Colv=FALSE, trace="none", dendrogram="none",
            hclustfun = .hclustWard,
            #col=blueyellow(100),
            #col=cm.colors(100),
            #col=hcl(h=240, c=30, l=seq(99,0,-1)),
            #col=hcl(h = seq(60, 258, by = 2)),
            #col=rgb.palette(100),
            col=lab.palette(100),
            #col=terrain.colors(100),
            #col=rainbow(100, end=.7),
            density.info=density, # "none", "density", "histogram"
            breaks=breaks,
            cexCol=.5)
  if (!is.null(fname)) dev.off()
}

.hclustWard <- function(x) {
  print(length(x))
  hclust(x, method="ward")
}

smoothImage <- function(data, margin, span=.05) {
  if (margin==1) {
    x <- 1:ncol(data)
  } else if (margin==2) {
    x <- 1:nrow(data)
  }
  out <- t(apply(data, margin, function(y) predict(loess(y~x, span=span))))
  return(out)
}

summarizeImageByCol <- function(data, col, chunkSize, FUN, span=NA, ...) {
  out <- foreach(i=isplitRows(data[,col], chunkSize=chunkSize), .combine="c") %do% {
    do.call(FUN, list(i, ...))
  }
  if (FUN=="bootCI") {
    out <- matrix(out, nrow=2, ncol=length(out)/2)
    if (!is.na(span)) {
      out <- apply(out, 1, function(x) predict(loess(x~c(1:length(x)), span=span)))
    }
  } else {
    if (!is.na(span)) {
      out <- predict(loess(out~c(1:length(out)), span=span))
    }
  }  
  return(out)
}

plotStatAndCI <- function(stat, CI, ...) {
  plot(stat, type="l", ...)
  a <- apply(CI, 2, lines, lty=2, ...)
}

orderByAnchor <- function(vals, anchor, width) {
  add_vector <- c(0, rep(c(1:width), each=2))
  add_vector <- add_vector * rep(c(1, -1), times=length(add_vector))
  ind_vector <- add_vector + anchor
 # return(ind_vector)
  order_vals <- order(vals[, ind_vector])
  return(order_vals)
  up <- c()
  
  #test <- c(vals[,1], vals[,2])
  #vals <- vals[order(-test),]
  return(vals[order_vals,])
                               
}

# Take 'image' matrix
# Smooth across columns
# Determine min, max
# Identify most 5' and 3' windows at some fraction of data range
identifyEnriched <- function(mat, span=0.05, fraction=.50) {
  mat <- smoothImage(mat, margin=1, span=span)
  mat <- t(apply(mat, 1, range01))
  #mat.max <- max(mat)
  #mat.min <- min(mat)
  #print(mat.max)
  #print(mat.min)
  #mat.range <- mat.max - mat.min
  #thresh <- mat.range * fraction
  #eps <- mat.range * .05
  #print(eps)
  #print(thresh)
  #mat <- t(apply(mat, 1, range01))
  mid <- round(ncol(mat) / 2)
  thresh_positions <- apply(mat, 1, function(x) {
    pos <- which(x>=(fraction-.05) & x<=(fraction+.05))    
    #pos <- which(x>=(thresh-eps) & x<=(thresh+eps)) 
    #if (is.null(length(pos))) {
    if (length(pos) == 0) {
      #return(list(min=(0 - mid), max=(0 + mid)))
      return(NA)
    }
    return(list(min=min(pos) - mid, max=max(pos) - mid))
  })
  #return(thresh_positions)
  return(thresh_positions[!is.na(thresh_positions)])
}

identifyEnrichedBed <- function(mat, span=0.05, fraction=.25, wsize=25, bed) {
  pos <- identifyEnriched(mat, span=span, fraction=fraction)
  #return(pos)
  pos <- lapply(pos, function(x) lapply(x, function(y) y * wsize))
  #return(pos)
  bed <- bed[match(names(pos), bed[,4]),]
  #return(bed)
  for (i in 1:nrow(bed)) {
    mid <- round((bed[i,3] + bed[i,2]) / 2)
    #print(mid)
    #print(pos[[i]])
    bed[i,2] <- mid + pos[[i]]$min
    bed[i,3] <- mid + pos[[i]]$max
  }    
  return(bed)
  
}
