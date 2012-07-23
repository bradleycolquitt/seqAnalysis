source("~/src/seqAnalysis/R/paths.R")
library(itertools)

## Generates matrix of samples values with annotations observations as rows and
##    positions as columns
##    Arguments:  data - profile data with [chr start stop name position strand value]
##                data_type - directory and normalization type
positionMatrix <- function(data, data_type="raw") {
  groups <- unique(data[,5])
  names <- unique(data[,4])
  data_groups = split(data, data[,5])
  out <- foreach(ind=icount(length(data_groups)), .combine="cbind") %dopar% {
    data_group <- data_groups[[ind]]
    m <- match(names, data_group[,4])
    return(data_group[m, 7])
  }
  rownames(out) <- names
  return(out)
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
      }
    }
  }  
}


makeImage <- function(sample, anno, data_type="unnorm", image=TRUE, heat=TRUE, range=NULL, average=NULL) {
  image.path <- paste(profile2.path, "norm", data_type, anno, "images", sample, sep="/")
  print(image.path)
  data <- as.matrix(read.delim(image.path, header=FALSE, row.names=1))
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
  return(data)
}


makeImage.all <- function(set="d3a", anno, data_type="unnorm/mean") {
  if (set=="d3a") {
    samples <- samples.d3a
  } else if (set=="cells") {
    samples <- samples.cells_norm
  }
  
  out <- lapply(samples, function(sample) {
    makeImage(sample, anno, data_type=data_type, image=FALSE)
  }) 
  names(out) <- samples
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

MP.heat <- function(data, density="none", range=NULL, average=NULL, fname=NULL) {
  require(gplots)
  if (is.null(fname)) {
    x11()
  } else {
    pdf(file=paste("~/s2/analysis/profiles/plots", fname, sep="/"), 12,12)
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
  rgb.palette <- colorRampPalette(c("yellow", "black", "blue"), space="rgb")
  lab.palette <- colorRampPalette(brewer.pal(5,"YlGnBu"), space="Lab")
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
            breaks=breaks)
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
  out <- apply(data, margin, function(y) predict(loess(y~x, span=span)))
  return(do.call("rbind", out))
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
MP.plyrMatrix <- function(mp) {
  fun <- function(x) {
    return(data.frame(x$ams_A, colnames=x$name))
  }
  out <- daply(mp, group, fun)
}
