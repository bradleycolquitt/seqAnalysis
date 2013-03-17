library(foreach)
library(plyr)

cuffdiff.path <- "~/s2/analysis/rna/cuffdiff"

CL.thresh <- function(de.file, FPKM=1, ratio=2, sig=TRUE, N=NULL, FDR=0) {
  de.data <- read.delim(de.file)
  
  ## Select signficants  
  if (sig & FDR == 0) {de.data <- subset(de.data, de.data$significant=="yes")
  } else if (!sig & FDR == 0) {de.data <- subset(de.data, de.data$significant=="no")
  } else {
    #q <- p.adjust(de.data$p_value, method="BH")
    #de.data$FDR <- q
    #de.data <- de.data[de.data$FDR <= FDR, ]
    de.data <- de.data[de.data$q_value <= FDR,]
  }
  ## Select genes with log2 FPKM above given value
  de.data <- de.data[log(de.data$value_1 + 1, 2) >= FPKM | log(de.data$value_2 + 1, 2) >= FPKM, ]
  
  ## Select genes with absolute ratio above given value
  #de.data <- transform(de.data, log2.FC = log(exp(ln.fold_change.), 2))
  de.data <- de.data[abs(de.data$log2.fold_change.) >= ratio, ]
  de.data$rel <- "down"
  de.data$rel[de.data$log2.fold_change. < 0] <- "up"

  #de.data$name2 <- id2names(de.data$test_id)
  de.data$strand <- id2strand(de.data$test_id)

  ## Combine isoforms
  de.data <- de.data[!duplicated(de.data$test_id),]
  
  if (!is.null(N)) {
    de.data.sp <- split(de.data, de.data$rel)
    tops <- lapply(de.data.sp, function(x) {
      q <- quantile(abs(x$log2.fold_change.), (nrow(x) - N)/nrow(x))
      x[abs(x$log2.fold_change.) >= q,]
    })
    de.data <- rbind(tops[[1]], tops[[2]]) 
  }
  
  return(de.data)
                                                   
}

CL.writeUpDown <- function(de.data, fname=NULL) {
  write.table(de.data[de.data$rel=="up",], paste(cuffdiff.path, paste(fname, "up", sep="_"), sep="/"),
              quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
  write.table(de.data[de.data$rel=="down",], paste(cuffdiff.path, paste(fname, "down", sep="_"), sep="/"),
              quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
}

CL.writeUpDownPositions <- function(de.data, fname=NULL) {
  up <- CL.getPositions(de.data[de.data$rel=="up",])
  down <- CL.getPositions(de.data[de.data$rel=="down",])
  write.table(up, paste(cuffdiff.path, paste(fname, "up", sep="_"), sep="/"),
              quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
  write.table(down, paste(cuffdiff.path, paste(fname, "down", sep="_"), sep="/"),
              quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
}

CL.getPositions <- function(de.data) {
  locus <- de.data$locus
  locus1 <- str_split(locus, ":")
  locus.chr <- unlist(lapply(locus1, function(x) x[1]))
  locus.pos <- unlist(lapply(locus1, function(x) x[2]))
  locus2 <- str_split(locus.pos, "-")
  locus.start <- unlist(lapply(locus2, function(x) x[1]))
  locus.end <- unlist(lapply(locus2, function(x) x[2]))
  locus.out <- data.frame(locus.chr,
                          as.numeric(locus.start),
                          as.numeric(locus.end),
                          as.character(de.data$test_id),
                          0,
                          as.character(de.data$strand))
  return(locus.out)
  
}

CL.volcano <- function(data, sample='value_1', filt.FDR=0.05, filt.ratio=1, title="", fname=NULL) {
  require(ggplot2)
  data <- data[log(data$value_1) >= 0 | log(data$value_2) >=0,]
  data.log2 <- log(data$value_2 / data$value_1,2)
  data.q <- -log(data$q_value)
  #data.q <- -log(p.adjust(data.p, method="BH"), 10)
  df <- data.frame(ratio=data.log2, pval=data.q)
  df <- df[is.finite(df$pval) & is.finite(df$ratio),]
  df$sig <- FALSE
  df$sig[df$pval >= filt.FDR & df$ratio >= filt.ratio |
         df$pval >= filt.FDR & df$ratio <= -filt.ratio] <- TRUE
#  return(df)
  if (is.null(fname)) {
    x11()
  } #else {
   #print(paste(plot.path, fname, sep="/"))
   # pdf(file=paste(plot.path, fname, sep="/"), 7, 7)
  #}  
  gp <- ggplot(df, aes(ratio, pval, color=sig))
  gp + geom_point(alpha=I(1/4)) + scale_color_manual("significant",
                    c("FALSE" = "black", "TRUE" = "red"))
  last_plot() + scale_x_continuous("log2 FPKM 2/FPKM1") + scale_y_continuous("-log10 q-values")
  print(last_plot() + opts(title=title,
                     panel.grid.minor = theme_blank(),
                     panel.grid.major = theme_blank(),
                     panel.background = theme_blank()))
  
}

RNA.boxplot <- function(data) {
  require(RColorBrewer)
  cols <- brewer.pal(3, "Dark2")
  ylim <- c(-10, 12)
#  layout(t(c(1,2)), widths=c(1,2))
                                        #par(mfrow=c(1,2), mar=c(2,1,2,0))
  bwx <- .4
  plot(1,1, xlim=c(0, 3), ylim=ylim, type="n", axes=FALSE, ann=FALSE)
  abline(h=0, col="light grey", lty=2)
  boxplot(data$value[!data$mk4], axes=FALSE, ylim=ylim, boxwex=bwx, col=cols[1], add=TRUE, at=1)
  axis(2)
  boxplot(value~hmc, data=data, subset=mk4 & hmc!="down", ylim=ylim, boxwex=.2, axes=F, at=c(1.5, 2.5), col=cols[2:3], add=TRUE)
  #abline(h=0, col="light grey", lty=2)
}

id2names <- function(id) {
  code <- read.delim("/seq/lib/roi/refgene_names")
  m <- match(id, code[,1])
  names <- code[m,2]
  return(names)
}

id2strand <- function(id) {
  code <- read.delim("/seq/lib/roi/refgene_names")
  m <- match(id, code[,2])
  strand <- code[m,3]
  return(strand)
}

RNA.filt <- function(data, thresh) {
  ind <- list(c(1,2), c(3,4), c(5,6))
  filt <- sapply(ind, function(i) data[,i[1]] >= thresh | data[,i[2]] >= thresh)
  data.filt <- data[as.logical(apply(filt, 1, prod)),]
  data.diff <- sapply(ind, function(i) data.filt[,i[2]] - data.filt[,i[1]])
  dimnames(data.diff) <- list(rownames(data.filt), c("5hmC", "5mC", "mRNA"))
  data.diff.r <- melt(data.diff)
  colnames(data.diff.r) <- c("Gene", "Measure", "Value")
  return(data.diff.r)
}


# Input:
#   1. matrix of enhancer-gene distances (rows - enhancers, columns - genes)
#   2. matrix of rna values

# Algorithm:
#   1. cut distances
#   2. for each enhancer, determine rna~cut distance ecdf for each rna group
#   3. average values for each rna group across enhancers


computeRnaEcdf <- function(eg_dist, rna, summary_stat="mean") {
#  eg_dist_cut <- t(apply(eg_dist, 1, cut, breaks=seq(0, by=20000, length.out=100)))
#  return(eg_dist_cut)
#  ecdf_all <- apply(eg_dist, 1, function(row) {
  ecdf_all <- foreach(row=isplitRows(eg_dist, chunkSize=1)) %do% {
    names(row) <- colnames(eg_dist)
    row_sorted <- sort(row)
    #return(names(row_sorted))
    row_cut <- cut(row_sorted, breaks=seq(0, by=2000, length.out=100))
    names(row_cut) <- names(row_sorted)
#    return(names(row_cut))
    row_ecdf <- apply(rna, 2, function(rna_col) {
      rna_col_select <- na.omit(rna_col[match(names(row_cut), rownames(rna))])
      return(rna_col_select)
#      return(cumsum(rna_col_select))
      #range01(cumsum(rna_col_select))
    })
    row_ecdf <- data.frame(cut=row_cut[match(rownames(row_ecdf), names(row_cut))], row_ecdf)
    return(row_ecdf)
  #})
  }
  ecdf_long <- do.call("rbind", ecdf_all)
#  return(ecdf_long)
  ecdf_summary <- ddply(ecdf_long, .(cut), function(d) apply(d[,2:ncol(d)], 2, summary_stat))
  ecdf_summary2 <- cbind(ecdf_summary[,1], apply(ecdf_summary[,2:4], 2, function(x) range01(cumsum(x))))
  return(ecdf_summary2)
}

# Data is 7 column bed: Column 1-6 are standard bed fields, Column 7 is measure
summarize_values_by_position <- function(data, window, step, FUN="mean") {
  chr_lengths <- read.delim("/seq/lib/mouse.mm9.genome", header=FALSE)
  data_split <- split(data, data[,1])
  values <- foreach(chrom=names(data_split), .combine="rbind") %dopar% {
    max_position <- chr_lengths[grep(paste(chrom, "$", sep=""), chr_lengths[,1]),2]
    data_curr <- data_split[[chrom]]
    values_curr <- foreach(pos=seq(1, max_position, step), .combine="c") %do% {
      subset <- data_curr[data_curr[,2]>=pos & data_curr[,3]<=(pos+window),]
      value <- do.call(FUN, list(subset[,7]))
      return(value)
    }
    return(data.frame(chrom=chrom, pos=1:length(values_curr), values=values_curr))
  }
  #names(values) <- names(data_split)
  return(values)
  values_df <- ldply(values, function(d) cbind(pos=1:length(d), value=d))
}
