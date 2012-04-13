

library(itertools)
library(foreach)
library(boot)
library(RColorBrewer)
library(multicore)
library(doMC)
library(MEDIPS)
library(cluster)
library(hopach)
library(gplots)
library(plyr)
library(reshape)
library(ggplot2)
library(proxy)
library(amap)
source("~/src/R/boot.R")
source("~/src/R/plotUtil.R")
source("~/src/seqAnalysis/R/paths.R")
source("~/src/MEDIPS/R/MEDIPS_mod.methylProfiling.R")
registerDoMC(cores=4)

#mean.cols <- c(brewer.pal(9, "Blues")[c(5,7,9)], brewer.pal(9, "Reds")[c(5,7,9)])

MP.feature <- function(data1, data2=NULL, feature=NULL, select=2, transf=FALSE, write=FALSE) {
  feature.data <- read.delim(feature,
                             header=FALSE)
  mp <- MEDIPS_mod.methylProfiling(data1=data1, data2=data2, ROI_file=feature.data, select=select, transf=transf)
  return(mp)
}


##
MP.feature.All <- function(data1, data2=NULL, set="general", select=2, transf=FALSE) {
  sample.name <- unlist(strsplit(data1@sample_name,".bed"))[1]
  mp.path.out <- paste(mp.path, sample.name, sep="/")
  if (!is.null(data2)){
    sample2.name <- unlist(strsplit(data2@sample_name, ".bed"))[1]
    mp.path.out <- paste(mp.path, paste(sample.name, sample2.name, sep="_"), sep="/")
  }
  if(!file.exists(mp.path.out)) dir.create(mp.path.out)
  if (set == "moe") feature.path <- feature.path.moe
  files.tbp <- checkExisting(feature.path, mp.path.out)
  print(files.tbp)
  x <- mclapply(files.tbp, mc.cores=3, mc.preschedule=FALSE, function(file) {
  #x <- lapply(files.tbp, function(file) {
    cat(paste("-- ", file, " --\n", sep=""))
    print(paste(mp.path.out, file, sep="/"))
    mp <- MP.feature(data1, data2, paste(feature.path, file, sep="/"), select=select, transf=transf)
    cat("Writing to file...\n")
    write.table(mp, file=paste(mp.path.out, file, sep="/"), quote=FALSE, sep="\t",
                row.names=TRUE, col.names=TRUE)  
    gc()
  })
  return(x)
}



MP.getIndVals.all <- function(files=NULL) {
  if (is.null(files)) files <- list.files(feature.path)
  vals <- foreach(file=files, .combine="rbind", .inorder=TRUE, .verbose=T) %dopar% {
    print(file)
    cells.hmc <- getIndVals(file, samples.hmc)
    cells.df.hmc <- makeDF(cells.hmc)
    cells.mc <-  getIndVals(file, samples.mc)
    cells.vals <- data.frame(ams=c(apply(cells.hmc, 2, rbind), apply(cells.mc, 2, rbind)),
                             cell=rep(cell.names,each=nrow(cells.hmc)),
                             DIP=rep(dip.names, each=nrow(cells.hmc) * 3),
                             feature=rep(file, times=nrow(cells.hmc) * 6))
    #cells.vals[,1:6] <- toNumeric(cells.vals[,1:6])
    return(cells.vals)
  }
  
  #vals[,1:6] <- apply(vals[,1:6], 2, function(x) as.numeric(x))
  return(vals)
}

## Extract feature regions by pval cut off (adjust first)
MP.thresh.p <- function(sample, feature, FDR=0.01, thresh=1, above=TRUE, write=FALSE) {
  data <- getFeatureFile(sample, feature)
  data <- data[is.finite(data$ratio),]
  data <- data[data$ratio > 0,]
  data.p <- data$p.value.wilcox
  data.q <- p.adjust(data.p, method="BH")
  data.thresh <- data[data.q <= FDR,]
  data.thresh <- data.thresh[data.thresh$rpm_A >= 0.15 | data.thresh$rpm_B >= 0.15,]
  #data.thresh <- data.thresh[data.thresh$rpm_A >= 0.1 | data.thresh$rpm_B >= 0.1,]
  #data.thresh <- data.thresh[data.thresh$ams_A >= 500 | data.thresh$ams_B >= 500,]
  if (above) {data.thresh <- data.thresh[log(data.thresh$ratio,2) >= thresh,]
  } else {data.thresh <- data.thresh[log(data.thresh$ratio,2) <= -thresh,]}

  if (write)  {
    lab <- NULL
    if (above) {lab <- "up"
    } else {lab <- "down"}
    out.path <- paste(mp.path, sample, paste(feature, "thresh", sep="_"), sep="/")
    if (!file.exists(out.path)) dir.create(out.path)
    write.table(data.thresh, file=paste(out.path,
                               paste(FDR, thresh, lab, sep="_"), sep="/"),
                quote=FALSE, sep="\t", row.names=TRUE, col.names=TRUE)
  }
  return (data.thresh)
}

MP.thresh.p.all <- function(dip, feature, FDR=0.01, thresh=1) {
  samples <- NULL
  if (dip == "hmedip") samples <- samples.hmc
  if (dip == "medip") samples <- samples.mc
  if (dip == "all") samples <- c(samples.hmc, samples.mc)
  if (dip == "d3a") samples <- samples.d3a
  foreach(sample= samples, .verbose=FALSE, .inorder=FALSE) %dopar% {
    features <- list.files(paste(mp.path, sample, sep="/"))
    features <- features[-grep("thresh", features)]
    #print(features)
    for(feature in features) {
      MP.thresh.p(sample, feature, FDR, thresh, above=TRUE, write=TRUE)
      MP.thresh.p(sample, feature, FDR, thresh, above=FALSE, write=TRUE)      
    }
  }  
}

MP.volcano <- function(sample, feature, title="", filt.FDR=2, filt.ratio=1, fname=NULL) {
  data <- getFeatureFile(sample, feature)
  data <- data[is.finite(data$ratio),]
  data.p <- data$p.value.wilcox
  data.q <- -log(p.adjust(data.p, method="BH"), 10)
  data.log2 <- log(data$ratio, 2)
  df <- data.frame(ratio=data.log2, pval=data.q)
  df <- df[is.finite(df$pval) & is.finite(df$ratio),]
  df$sig <- FALSE
  df$sig[df$pval >= filt.FDR & df$ratio >= filt.ratio |
         df$pval >= filt.FDR & df$ratio <= -filt.ratio] <- TRUE
  if (is.null(fname)) {
    x11()
  } else {
    print(paste(plot.path, fname, sep="/"))
    pdf(file=paste(plot.path, fname, sep="/"), 7, 7)
  }  
  gp <- ggplot(df, aes(ratio, pval, color=sig))
  gp + geom_point(alpha=I(1/4)) + scale_color_manual("sig",
                    c("FALSE" = "black", "TRUE" = "red"))
  last_plot() + scale_x_continuous("log2 5hmC") + scale_y_continuous("-log10 adjusted p-values")
  print(last_plot() + opts(title=title,
                     panel.grid.minor = theme_blank(),
                     panel.grid.major = theme_blank(),
                     panel.background = theme_blank()))
  
  if (!is.null(fname)) {dev.off()}
                                        #sub_data <- subset(df, df$ratio>=2 | df$ratio<=2 & df$pval >= 2)
  #gp_sub <- gpplot()
                                        #points(sub_data$ratio,sub_data$pval, col="red", pch=20)
  #last_plot() + geom_point(aes(ratio[ratio>=2 | ratio<=-2 &], pval[pval>=2], color="red"))


  #qplot(ratio, pval, data=df, alpha=I(1/5))
  #plot(data.log2, data.q, pch=20)
}

MP.makeFeatureSummaries <- function(feature, sample, sample.name) {
  out.path <- "~/storage/analysis/mprofiles/features/feature_summaries"
  #sample.name <- NULL
  #if (sample == samples.mc) sample.name <- "mc"
  vals <- getIndVals(feature, sample)
  write.table(vals, file=paste(out.path, paste(sample.name, feature, sep="_"), sep="/"),
              quote=FALSE, sep="\t")  
}

MP.makeFeatureSummaries.2 <- function(feature, sample_group = "d3a") {
  if (sample_group == "d3a") {
    data.mc <- getTwoVals(feature, "moe_wt_mc_moe_dnmt3a_mc")
    data.hmc <- getTwoVals(feature, "moe_wt_hmc_moe_dnmt3a_hmc")
  }
  out <- rbind(data.mc, data.hmc)
  #return(out)
  out <- data.frame(Feature=rep(feature, times=nrow(out)), Id=rownames(out), Modification=rep(c("5mC", "5hmC"), each=nrow(data.mc)), WT=out[,1], KO=out[,2])
  out <- melt(out, measure=c("WT", "KO"))
  return(out) 
}

MP.makeFeatureSummaries.All <- function(sample_group = "d3a") {
  out.path <- "~/storage/analysis/mprofiles/features/feature_summaries"
  features <- list.files(feature.path)
  data <- foreach(feature=features, .combine="rbind", .verbose=TRUE) %dopar% {
    return(MP.makeFeatureSummaries.2(feature, sample_group))
  }
  return(data)
}

## Construct matrix of given value from given features
MP.makeAMSmatrix <- function(feature, samples="all", N=200, value_type="rpm") {
  vals <- NULL
  if (samples=="all") {
    hmc.val <- getIndVals(feature, samples.hmc, value_type)
    mc.val <- getIndVals(feature, samples.mc, value_type)
    vals <- cbind(hmc.val, mc.val)
  } else if (samples == "hmc") {
    vals <- getIndVals(feature, samples.hmc)
  } else if (samples == "mc") {
    vals <- getIndVals(feature, samples.mc)
  } else if (samples == "sc") {
    vals <- getIndVals(feature, samples.sc)
  } else if (samples == "ac3") {
    tmp.val <- getFeatureFile(samples.ac3, feature)
    vals <- cbind(tmp.val$ams_A, tmp.val$ams_B)
    colnames(vals) <- c("WT", "AC3")
  } else if (samples == "d3a") {
    vals.hmc <- getTwoVals(feature, samples.d3a[1], value_type)
    vals.mc <- getTwoVals(feature, samples.d3a[2], value_type)
    vals <- cbind(vals.hmc, vals.mc)
    colnames(vals) <- c("5hmC WT", "5hmC KO", "5mC WT", "5mC KO")
  }
  #return(vals)
  vals <- threshRowsByVal(vals, 300)
  cells.var <- unlist(apply(vals, 1, var))
  index <- c(1:nrow(vals))
  if (N > 0) index <- MP.thresh.byQ(cells.var, N=N)
  return(vals[index,])
}

MP.heatmap.2 <- function(vals, rowv=NULL, trim_dist=FALSE, lab=NULL, fname=NULL, ...) {
  if (is.null(fname)) {X11("", 10, 10)
  } else {pdf(file=paste(cluster.path, "heatmap", paste(fname, ".pdf", sep=""), sep="/"),
              width=12, height=12)}             
  if (trim_dist) {
    distfun <- distfun_trim
  } else {
    distfun <- function(x) {dist(x, method="Euclidean")}
    #distfun <- dist
  }
  if (!is.null(rowv)) {
    print("rowv")
    rowv <- rowv
  } else {
    rowv <- TRUE
  }

  heatmap.2(vals, Rowv=rowv, Colv=FALSE, dendrogram="row", trace="none",
            #col=topo.colors(100),
            distfun=distfun,
            hclustfun=.hclustWard,
            col=greenred(100),
            #labCol=colnames(vals),
            labCol=lab,
            labRow="",
            density.info="density",
            keysize=.75,
            margins=c(5,5), ...)
  if(!is.null(fname)) {dev.off()}        
}

distfun_trim <- function(x) {
  print(dim(x))
  dist(x[,c(1:3)])
}

.hclustWard <- function(x) {
  print(length(x))
  hclust(x, method="ward")
}

MP.makeSigMatrix <- function(feature, thresh=700) {
  data.ex <- getFeatureFile(samples[1], feature)
  data.mat <- matrix(0, nrow=nrow(data.ex), ncol=length(samples) * 2)
  rownames(data.mat) <- rownames(data.ex)
  for(i in 1:length(samples)) {
    data.up <- MP.thresh.p(samples[i], feature, thresh=thresh, above=TRUE)
    data.mat[rownames(data.mat)%in%rownames(data.up),i] <- 1
    data.down <- MP.thresh.p(samples[i], feature, thresh=thresh, above=FALSE)
    data.mat[rownames(data.mat)%in%rownames(data.down),i+1] <- 1
  }
  data.mat <-data.mat[rowSum(data.mat) > 0,]
  return(data.mat)
}
MP.clusterAMS.PAM <- function(vals) {
  result <- list()
  for(i in 2:7) {
    result <- c(result, list(pam(vals, k=i)))
  }
  return(result)
}
MP.cluster.clara <- function(vals) {
  result <- list()
  for(i in 2:7) {
    result <- c(result, list(clara(vals, k=i)))
  }
  return(result)
}
MP.plot.PAM <- function(pams, fname=NULL) {
  if (is.null(fname)) {x11("", 10, 10)
  } else {png(file=fname)}             
  par(mfrow=c(ceiling(length(pams)/3), 3))
  x <- lapply(pams, function(x) plot(x, which.plot=1))
  if (!is.null(fname)) dev.off()
}
MP.PAM.saveByCluster <- function(pam, fname) {
  cluster.path <- "~/analysis/mprofiles/features/clustering/pam"
  file.path <- paste(cluster.path, fname, sep="/")
  clusters <- pam$clustering
  for (i in 1:nrow(pam$clusinfo)) {
    write.table(names(clusters)[clusters==i], file=paste(file.path, i, sep="_"),
                quote=FALSE, row.names=FALSE, col.names=FALSE)
  }
}
MP.clusterAMS.hclust <- function(vals, groups, fname=NULL) {
  vals.dist <- dist(vals)
  result <- hclust(vals.dist, method="ward")
  if (groups == 0) return(result)
  result.cut <- cutree(result, groups)
  if(!is.null(fname)) {
    for(i in 1:groups) {
      ind <- rownames(vals) %in% names(result.cut)[result.cut==i] 
      out <- vals[ind,]
      write.table(out,
                  file=paste(cluster.path, "hclust",
                  paste(fname, "cut", groups, i, sep="_"), sep="/"),
                  quote=FALSE, row.names=TRUE, col.names=TRUE, sep="\t")
    }
  }
  return(result.cut)
  
}

MP.clusterAMS.hopach <- function(vals) {
  vals.dist <- distancematrix(vals, "cosangle")
  vals.obj <- hopach(vals, dmat=vals.dist)
  return(list(vals.dist,vals.obj))
}

