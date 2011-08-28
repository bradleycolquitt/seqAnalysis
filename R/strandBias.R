classes <- c("character", "numeric", "numeric", "character", "character", "numeric")
strand.path <- "~/storage/analysis/strand"
samples <- paste(rep(c("omp", "ngn", "icam"), each=2), c("hmedip", "medip"), sep="_")

strandBias <- function(sample, roi) {
  plus <- paste(strand.path, sample, paste(roi, "plus", sep="_"), sep="/")
  minus <- paste(strand.path, sample, paste(roi, "minus", sep="_"), sep="/")
  plus.data <- read.delim(plus, header=FALSE, colClasses=classes)
  minus.data <- read.delim(minus, header=FALSE, colClasses=classes)
  scores <- vector(length=nrow(plus.data))
  plus.ind <- plus.data[,4] == "+"
  minus.ind <- plus.data[,4] == "-"
  totals <- plus.data[,6] + minus.data[,6]
  scores[plus.ind] <- (plus.data[plus.ind,6] - minus.data[plus.ind,6])/totals[plus.ind]
  scores[minus.ind] <- (minus.data[minus.ind,6] - plus.data[minus.ind,6])/totals[minus.ind]

  out.data <- data.frame(chr=plus.data[,1],
                         start=plus.data[,2],
                         end=plus.data[,3],
                         strand=plus.data[,4],
                         name=plus.data[,5],
                         plus=plus.data[,6],
                         minus=minus.data[,6],
                         total=totals,
                         score=scores)

  write.table(out.data, file=paste(strand.path, sample,
                          paste(roi, "analyzed", sep="_"), sep="/"),
              quote=FALSE, sep="\t", row.names=FALSE)
  return(out.data)
}

SB.AllSamples <- function(roi) {
  for (sample in samples) {
    x <- strandBias(sample, roi)
  }
}

SB.summarize <- function(roi) {
  samp <- list()
  for (sample in samples) {
    data <- read.delim(paste(strand.path, sample, paste(roi, "analyzed", sep="_"), sep="/"), as.is=TRUE)
    samp <- c(samp, sample=list(data)) 
  }
  names(samp) <- samples
  
  chr <- SB.reform(samp, 'chr')
  strand <- SB.reform(samp, 'strand')
  name <- SB.reform(samp, 'name')
  total <- SB.reform(samp, 'total')
  scores <- SB.reform(samp, 'score')
  samples.ext <- rep(names(samp), each=nrow(samp[[1]]))
  samples.split <- str_split(samples.ext, "_")
  cells <- unlist(lapply(samples.split, function(x) x[1]))
  dips <- unlist(lapply(samples.split, function(x) x[2]))
  
  rna <- read.delim("~/storage/data/rna/cuffdiff/refgene_name2/genes.fpkm_tracking")                                        
  rna.ind <- match(samp[[1]]$name, rna[,1])
  rna.omp <- rna[rna.ind, 'q1_FPKM']
  rna.ngn <- rna[rna.ind, 'q2_FPKM']
  rna.val <- c(rep(rna.omp, times=6), rep(rna.ngn, times=6))
  or <- rep(FALSE, times=length(name))
  or[grep("Olfr", name)] <- TRUE
            
  out <- data.frame(chr=as.character(chr),
                    strand=as.character(strand),
                    name=as.character(name),
                    total=total,
                    score=scores,
                    cell=factor(cells),
                    dip=factor(dips),
                    or=factor(or),
                    rna.sample=rep(c("omp", "ngn"), each=6*nrow(samp[[1]])),
                    rna.log2=log(rna.val, 2))
  write.table(out, file=paste(strand.path, "df",
                          paste(roi, "summary", sep="_"), sep="/"),
              quote=FALSE, sep="\t", row.names=FALSE)
  return(out)
}

strand_path <- "~/s2/analysis/strand"
samples_d3a <- paste("moe", rep(c("wt", "d3a"), each=3), rep(c("mc", "hmc", "in"), each=2), sep="_")
samples_cells <- paste(rep(c("omp", "ngn", "icam"), each=2), c("medip", "hmedip"), sep="_")
SB.getData <- function(sample_group, roi, group=FALSE) {
  require(reshape)
  if (sample_group == "cells") {
    samples <- samples_cells
  } else if(sample_group == "d3a") {
    samples <-  samples_d3a
  }
  data <- lapply(samples, function(x) {
    path <-paste(strand_path, x, paste(roi, "_scores.txt", sep=""), sep="/")
    return(read.delim(path, header=FALSE))
                      #colClasses=c('character', 'numeric', 'numeric', 'character', 'numeric', 'character', 'numeric')))
  })
  names(data) <- samples
  data.all <- ldply(data, data.frame)
  colnames(data.all) <- c("sample", "chr", "start", "end", "name", "blank", "strand", "total", "score")
  lib <- unique(data.all$sample)
  split_lib <- str_split(lib, "_")
    if (sample_group == "cells") {
    label_ind <-  c(1,2)
    } else if (sample_group == "d3a") {
    label_ind <- c(2,3)
    }
  matchLib <- function(lib, split_lib, labs) {
    out_labs <- split_lib[match(labs, lib)] 
    return(out_labs)
  }
  data.all <- cbind(sample = unlist(lapply(matchLib(lib, split_lib, data.all$sample), function(x) x[label_ind[1]])), ip = unlist(lapply(matchLib(lib, split_lib, data.all$sample), function(x) x[label_ind[2]])), data.all[,2:ncol(data.all)])
  if (!group) {
    return(data.all)
  } else {
    df <- data.frame(nrow=nrow(data[[1]]), ncol=2+length(data))
                                        #return(df)
    df <- data.frame(data[[1]][,4], data[[1]][,5])
    for (i in 1:length(data)) {
      df[,i+2] <- data[[i]][,7]
    }
    colnames(df) <- c("id", "pos", names(data))
    return(df)
  }  
}

SB.makeDensities <- function(data, group=FALSE) {
  require(ggplot2)
  x11()
  gg <- ggplot(data, aes(score, color=sample))
  gg + geom_density() + facet_grid(.~ip)
}
generateProfiles <- function(samples, roi) {
  data <- SB.getData(samples, roi)
  data.mean <- lapply(data, function(x) {
    return(tapply(x[,7], x[,5], mean))
  })
  
  return(data.mean)
}


SB.computeBD <- function(data) {
  vals <- apply(data, 1, maxMin)
  return(vals)
}


maxMin <- function(vec) {
  return(max(vec) - min(vec))
}

SB.reform <- function(data.list, select) {
  x <- foreach(i=data.list, .combine="c") %do% {
    return(i[,select])
  }
  return(x)
}
