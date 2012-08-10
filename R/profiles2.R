
library(itertools)
library(foreach)
library(boot)
library(RColorBrewer)
library(multicore)
library(doMC)
library(MEDIPS)
library(colorspace)
library(cluster)
library(gplots)

source("~/src/seqAnalysis/R/paths.R")
source("~/src/seqAnalysis/R/boot.R")
source("~/src/seqAnalysis/R/plotUtil.R")
source("~/src/MEDIPS/R/MEDIPS_mod.methylProfiling.R")

registerDoMC(cores=4)

col4 <- brewer.pal(4, "Spectral")
col3 <- brewer.pal(3, "Dark2")
col2 <- rev(brewer.pal(3, "Set1")[1:2])

# Accepts data from makeProfile and writes to givin out.path
writeModule <- function(out.path, sample, group2=NULL, fun, rm.outliers=0, profile, CI) {
  if (!file.exists(out.path)) {
    dir.create(out.path)
  }
  out.name <- NULL

  if (is.null(group2)) {
    out.name <- sample
  } else {
    out.name <- paste(sample, group2, sep="_")
  }
  if (rm.outliers == 0) {
    out.name <- out.name
  } else {
    out.name <- paste(out.name, paste("trim", rm.outliers, sep=""), sep="_")
  }
  
  # Test if data list has multiple elements, indicating that profile has been split
  if (length(profile) > 1) {
    lapply(names(profile), function(name) {
      write.table(as.matrix(profile[[name]]),
                  #file=paste(out.path, paste(out.name, name, "mean", sep="_"), sep="/"),
                  file=paste(out.path, paste(out.name, name, fun, sep="_"), sep="/"),
                  quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
      lapply(c(1:2), function(x) {
        unwrapped <- NULL
        if (is.null(group2)) {unwrapped <- .unwrap(CI[[name]], x)
        } else {unwrapped <- apply(CI[[name]], 2, function(y) .unwrap(y, x))}
        write.table(unwrapped,
                  #file=paste(out.path, paste(out.name, name, "mean_bootCI", x, sep="_"),
                    file=paste(out.path, paste(out.name, name, fun, "bootCI", x, sep="_"),
                    sep="/"),                    
                    quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
      })
    })
  } else {
    write.table(profile,
              #file=paste(out.path, paste(out.name, "mean", sep="_"), sep="/"),
                file=paste(out.path, paste(out.name, fun, sep="_"), sep="/"),
              quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)     
    
    lapply(c(1:2), function(x) {
      unwrapped <- NULL
      if (is.null(group2)) {unwrapped <- .unwrap(CI[[1]], x)
      } else {unwrapped <- apply(CI[[1]], 2, function(y) .unwrap(y, x))}
      write.table(unwrapped,
                    #file=paste(out.path, paste(out.name, "mean_bootCI", x, sep="_"),
                    file=paste(out.path, paste(out.name, fun, "bootCI", x, sep="_"),
                    sep="/"),
                    quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
    })
  }  
}


profileCompute <- function(data, param) {
  grouping <- lapply(param[[2]], function(x) if(!is.null(x)) data[,x])
  if (length(param) == 3) {
    return(tapply(data[,param[[1]]], grouping, param[[3]]))
  } else {
    return(tapply(data[,param[[1]]], grouping, param[[3]], param[[4]]))
  }  
}

## Takes intersection information produced by group_2.py and generates meta-plot of sample
##    over given annotation.
## Arguments: anno - annotation file (from ~/lib/)
##            sample - sequencing sample whose read per window values were interseted with anno
##            group2 - optional, group rows of anno by some factor to process groups independently
##            data_type - type of normalization. also directory where data is located
##            rm.outliers - fraction of extreme-valued samples to remove from each position
##                          e.g. 0.05 removes top and bottom 5% of values
##            sample_num - number of annotation observations to sample. 0 indicates to use whole set
##            write - write data to file
makeProfile2 <- function(anno, samples, group2=NULL, data_type="unnorm/mean", fun="mean", rm.outliers=0, sample_num=0, write=TRUE) {
  sample_path <- paste(profile2.path, "norm", data_type, anno, sep="/")
  a <- sapply(samples, function(sample) {
      data <- read.delim(paste(sample_path, sample, sep="/"))
      colnames(data) <- c("chr", "start", "end", "name", "group", "strand", "norm")
      ind <- "norm"
      profile <- NULL
      CI <- NULL
      if (sample_num > 0) {
        print(sample_num)
        sample_names <- sample(unique(data$name), sample_num)
        data <- data[data$name %in% sample_names,]
        gc()
      }

      if (rm.outliers > 0) {
        thresh <- quantile(data$norm, probs=c(rm.outliers, 1-rm.outliers))
        data <- data[data$norm >= thresh[1] & data$norm <= thresh[2], ]
      }

      if (!is.null(group2)) {
        group2_vals <- read.delim(paste(group2.path, group2, sep="/"), header=FALSE)
        data$group2 <- group2_vals[match(data$name, group2_vals[,1]),2]
        group2_name <- "group2"
        #profile <- lapply(ind, function(x) profileCompute(data, list(x, list("group", group2_name), mean)))
        profile <- lapply(ind, function(x) profileCompute(data, list(x, list("group", group2_name), fun)))
        names(profile) <- ind
        #CI <- lapply(ind, function(x) profileCompute(data, list(x, list("group", group2_name), bootCI)))
        CI <- lapply(ind, function(x) profileCompute(data, list(x, list("group", group2_name), bootCI, fun)))
        names(CI) <- ind
      } else {
        #profile <- lapply(ind, function(x) profileCompute(data, list(x, list("group"), mean)))
        profile <- lapply(ind, function(x) profileCompute(data, list(x, list("group"), fun)))
        names(profile) <- ind
        #CI <- lapply(ind, function(x) profileCompute(data, list(x, list("group"), bootCI)))
        CI <- lapply(ind, function(x) profileCompute(data, list(x, list("group"), bootCI, fun)))
        names(CI) <- ind
      }
      #print("Finished compute")
      if(!write) {
        return(profile)
      } else {
        #out.path <- paste(anno, "profiles", sep="/")
        out.path <- paste(sample_path, "profiles", sep="/")
        group2_out <- NULL
        if (!is.null(group2)) {
          group2_out <- str_split(group2, "/")
          group2_out <- group2_out[length(group2_out)]
        }
        
        writeModule(out.path=out.path, sample=sample, group2=group2_out, fun=fun,
                    rm.outliers=rm.outliers, profile=profile, CI=CI)   
      }
      rm(profile)
      gc()
  })
 
}

## Wrapper to run makeProfile2 on all samples within a given profile direction
## Arguments as makeProfile2
makeProfile2.allSamp <- function(anno, group2=NULL, data_type="unnorm/mean", fun="mean", rm.outliers=0, sample_num=0, write=T) {  
  sample_path <- paste(profile2.path, "norm", data_type, anno, sep="/")
  print(sample_path)
  samples <- list.files(sample_path)
  ind <- grep("profiles", samples)
  if (length(ind) > 0) samples <- samples[-ind]
  out_path = paste(sample_path, "profiles", sep="/")
  data <- foreach(sample=samples) %do% {
    out_name <- sample
    if (!is.null(group2)) {
      out_name <- paste(out_name, group2, sep="_")
    }
    if (rm.outliers > 0) {
      out_name <- paste(out_name, paste("trim", rm.outliers, sep=""), sep="_")
    }
    if (file.exists(paste(out_path, paste(out_name, "mean", sep="_"), sep="/"))) {
      print("Skipping")
      next 
    }
    print(sample)
    return(makeProfile2(anno, sample, data_type=data_type, fun=fun, group2=group2,
                           rm.outliers=rm.outliers, sample_num=sample_num, write=write))
  }
}

## Wrapper to run makeProfile2 on all annotation and all samples within given profile data_type
## Arguments as makeProfile2
makeProfile2.allAnno <- function(data_type="rpm_avg_2", group2=NULL, write=TRUE) {
  files <- list.files(paste(profile2.path, "norm", data_type, sep="/"))
  registerDoMC(cores=4)
  i <- 0
  for(file in files) {
    i <- i + 1
    cat("--")
    cat(file)
    cat(paste(" ", i, " of ", length(files), sep=""))
    cat("\n")
    tryCatch(makeProfile2.allSamp(file, group2=group2, data_type=data_type, write=write), error = function(e) {
      print(paste("Skipping", file, sep=" "))
      print(e)
      return
    })
  } 
}

## Draws profile and confidence intervals
plotProfiles <- function(profile=NULL, ci_1=NULL, ci_2=NULL, cols=NULL,
                            smooth=FALSE, span=0.1) {

  col.rgb <- col2rgb(cols, alpha=TRUE)
  col.rgb[4,] <- 50
  
  if (smooth) {
    profile <- smoothProfile(profile, span)
    if (length(ci_1) > 0) {
      ci_1 <- smoothCI(ci_1, span)
      ci_2 <- smoothCI(ci_2, span)
    }
  }

  # if CI provided, draw polygon
  if (!is.null(ci_1)) {
      polygon(c(c(1:length(profile)), c(length(profile):1)),
              c(ci_1, rev(ci_2)),
              col=rgb(col.rgb[1,1], col.rgb[2,1], col.rgb[3,1],
                alpha=col.rgb[4,1], maxColorValue=255),
              border=F)
  }
  lines(profile, col=cols)
}

## Establishes plot area
plotAnno <- function(data, annotation, wsize, cols=NULL, lab=NULL,
                        y.val=NULL, combine=FALSE,
                        stack=FALSE, ...) {
  lab.data <- NULL
  if (is.null(dim(data[[1]][[1]]))) {
    x.lim <- length(data[[1]][[1]])
  } else {
    x.lim <- nrow(data[[1]][[1]])
  }
 
  
  #if (stack || length(data) > 1) {
  # Set up plot area
  plot(1, 1, type="n",
       xlim=c(1, x.lim),
       #ylim=c(round(y.val[1],2), round(y.val[2], 2)),
       ylim=y.val,
       xlab="",
       ylab="",
       ann=FALSE, axes=FALSE)

  # Send data to drawer  
  for(i in 1:length(data)) {
      data.curr <- data[[i]]
      data.val <- data.curr[[1]]
      data.ci1 <- data.curr[[2]]
      data.ci2 <- data.curr[[3]]
      if (ncol(data.val) > 1) {
        for (i in 1:ncol(data.val)) {
          plotProfiles(profile=data.val[,i], ci_1=data.ci1[,i], ci_2=data.ci2[,i],
                          smooth=FALSE, cols=cols[i])  
        }
        lab.data <- computeAxis(data.val[,1], wsize, lab)
      } else {
        plotProfiles(profile=data.val[,1], ci_1=data.ci1[,1], ci_2=data.ci2[,1],
                          smooth=FALSE, cols=cols[i])
        lab.data <- computeAxis(data.val[,1], wsize, lab)
      }
    }
  
  ## Draw and label axes
  box()
  axis(1, at=lab.data$pos,
       labels=c(paste("-", lab.data$dist, " kb", sep=""), lab,
                      paste("+", lab.data$dist, " kb", sep="")),
       cex.axis=1)
  print(y.val)
  axis(2, at=y.val)
  #axis(2, at=seq(round(y.val[1], 1), round(y.val[2],1), round(diff(y.val)/3, 2)), cex.axis=1)
  abline(v=lab.data$pos[2], lty=2, col="grey")
  if (length(lab) == 2) abline(v=lab.data$pos[3], lty=2, col="grey")
}

##############################################################
## Primary interface for drawing profile of a single sample ##
##############################################################
plot2 <- function(annotation, sample, orient=2, data_type="unnorm/mean", fun="mean", group2=NULL,
                     cols=1, lab=c("",""), y.vals=NULL, wsize=25, range=NULL, type="range", fname=NULL) {

  ## make profile path
  if (!is.null(group2)) {
    anno <- paste(annotation, group2, data_type, sep="_")
  } else {
    anno <-  paste(annotation, data_type, sep="_")
  }
  
  ## Orient 1 data structure is no longer used (2/14/12)
  if (orient==1) {
    data <- profileRead(paste(profile2.path, sample, "profiles", sep="/"), fun,
                        paste(annotation, group2, data_type, sep="_"))
  } else if (orient==2) {
    data <- profileRead(paste(profile2.path, "norm", data_type, annotation, "profiles", sep="/"), fun, sample, group2)
  }

  ## If range specificed, trim profile by given values
  if (!is.null(range)) {
      data <- lapply(data, function(x) x[range[1]:range[2],])
  }

  ## Set up graphics device
  if (is.null(fname))  {
    x11("", 6, 4)
  } else if (fname=="manual") {
    # do nothing, device has already been set
  } else {
    if (fname=="auto")
      dt <- paste(unlist(str_split(data_type, "/")), collapse="_")
      if (!is.null(group2)) {
        fname <- paste(annotation, sample, dt, group2, sep="_")
      } else {
        fname <- paste(annotation, sample, dt, sep="_")        
      }
    fname <- paste(fname, ".pdf", sep="")
    print(paste("Saving to ", fname, sep=""))
    pdf(file=paste(profile2.path, "norm", "plots", fname, sep="/"), 6, 4.5)
  }
  
  ## If y.vals not provided, determine y axis range from the data
  if (is.null(y.vals)) y.vals <- getRange(list(data), buffer=0)
  
  ## Send data to plotAnno
  plotAnno(list(data), annotation, cols=cols, lab=lab,
                y.val=y.vals, wsize=wsize) 
  
  ## Label axes
  #mtext("AMS", at=.775, side=2, outer=T, line=-1, cex=1.6)
  #mtext("AMS", at=.275, side=2, outer=T, line=-1, cex=1.6)
  #mtext("Normalized read count", at=.5, side=2, outer=T, line=-2, cex=1)
  #par(las=1)
  #mtext("5hmC", at=.775, side=2, outer=T, line=1, cex=1.6)
  #mtext("5mC", at=.275, side=2, outer=T, line=1, cex=1.6)
  #mtext("5hmC", at=.275, side=3, outer=T, line=0, cex=1.6)
  #mtext("5mC", at=.775, side=3, outer=T, line=0, cex=1.6)
  #mtext(paste(sample, " -> ", group2, sep=""), side=3, outer=T, cex=1.6)
  if (!is.null(fname)) {
    if (fname != "manual") dev.off()
  }
}

plot2.several <- function(annotation, set="d3a", data_type="unnorm/mean", group2=NULL, cols=NULL, lab=c("",""), y.vals=NULL, wsize=25, standard=FALSE, range=NULL, baseline=FALSE, fname=NULL) {
  samples <- NULL
  orient <- 1
  rows <- 1
  columns <- 1
  if (set=="d3a") {
    samples <- list(list("moe_wt_hmc.bed", "moe_d3a_hmc.bed"), list("moe_wt_mc.bed", "moe_d3a_mc.bed"))
    legend <- c("Dnmt3a +/+", "Dnmt3a -/-")
    rows <- 2
    columns <- 1
  } else if (set=="d3a_2") {
    samples <- list(list("moe_d3a_wt_hmc", "moe_d3a_ko_hmc"), list("moe_d3a_wt_mc", "moe_d3a_ko_mc"))
    rows <- 2
    columns <- 1
    orient <- 2
  } else if (set=="d3a_3") {
    samples <- list(list("d3a_wt_hmc", "d3a_ko_hmc"), list("d3a_wt_mc", "d3a_ko_mc"))
    rows <- 2
    columns <- 1
    orient <- 2
  } else if (set=="d3a_4") {
    samples <- list(list("moe_d3a_wt_hmc_30M_rpkm",  "moe_d3a_ko_hmc_rpkm"),
                    list("moe_d3a_wt_mc_rpkm", "moe_d3a_ko_mc_rpkm"))
    rows <- 2
    columns <- 1
    orient <- 2
  } else if (set=="d3a_rlm") {
    samples <- list(list("moe_wt_hmc_rlm", "moe_d3a_hmc_rlm"), list("moe_wt_mc_rlm", "moe_d3a_mc_rlm"))
    rows <- 2
    columns <- 1
  } else if (set=="d3a_rf") {
    samples <- list(list("moe_wt_mc_rf", "moe_d3a_mc_rlm"))
    rows <- 1
    columns <- 1
  }
  else if (set=="cells") {
    samples <- list(list("omp_hmedip.bed", "ngn_hmedip.bed", "icam_hmedip.bed"),
                    list("omp_medip.bed", "ngn_medip.bed", "icam_medip.bed"))
    rows <- 2
    columns <- 1
    legend <- c("mOSN", "GBC", "HBC")
  } else if (set=="cells_rlm") {
    samples <- list(list("omp_hmc_rlm", "ngn_hmc_rlm", "icam_hmc_rlm"),
                    list("omp_mc_rlm", "ngn_mc_rlm", "icam_mc_rlm"))
    rows <- 2
    columns <- 1
  } else if (set=="cells_norm") {
    samples <- list(list("omp_hmc", "ngn_hmc", "icam_hmc"),
                    list("omp_mc", "ngn_mc", "icam_mc"))
    columns <- 1
    #if (!is.null(group2)) {
    #  columns <- 3
    #}
    rows <- 2
    orient <- 2
  } else if (set=="cells_rpkm") {
    samples <- list(list("omp_hmc_rpkm", "ngn_hmc_rpkm", "icam_hmc_rpkm"),
                    list("omp_mc_rpkm", "ngn_mc_rpkm", "icam_mc_rpkm"))
    columns <- 1
    #if (!is.null(group2)) {
    #  columns <- 3
    #}
    rows <- 2
    orient <- 2
  } else if (set=="cells_norm_hmc") {
    samples <- list("omp_hmc", "ngn_hmc", "icam_hmc")
    columns <- 1
    rows <- 2
    orient <- 2
  } else if (set=="cells_norm_mc") {
    samples <- list("omp_mc", "ngn_mc", "icam_mc")
    columns <- 1
    rows <- 2
    orient <- 2
  } else if (set=="tfo_omp") {
    samples <- list(list("tfo_hmc", "omp_hmc"),
                    list("tfo_mc", "omp_mc"))
    rows <- 2
    columns <- 1
    orient <- 2
  } else if (set=="tfo_omp_rpkm") {
    samples <- list(list("tfo_hmc_22M", "omp_hmc_rpkm"))
    rows <- 1
    columns <- 1
    orient <- 2
  } else if (set=="tfo_omp_d3a") {
    samples <- list(list("tfo_hmc", "omp_hmc", "moe_d3a_wt_hmc", "moe_d3a_ko_hmc"),
                    list("tfo_mc", "omp_mc", "moe_d3a_wt_mc", "moe_d3a_ko_mc"))
    rows <- 2
    columns <- 1
    orient <- 2
  } else if (set=="nuc") {
    samples <- list(list("omp_nuc_0123", "icam_nuc_01234"))
    rows <- 1
    columns <- 1
    orient <- 2
  } else if (set=="tt3") {
    samples <- list(list("o.tt3.1_hmc_rpkm", "o.tt3.2_hmc_rpkm"),
                    list("o.tt3.1_mc_rpkm", "o.tt3.2_mc_rpkm"))
    rows <- 2
    columns <- 1
    orient <- 2
  } else if (set=="tt3_2") {
    samples <- list(list("omp_hmc_rpkm", "o.tt3.2_hmc_rpkm"),
                    list("omp_mc_rpkm", "o.tt3.2_mc_rpkm"))
    rows <- 2
    columns <- 1
    orient <- 2
  } else if (set=="d3a_nuc") {
    samples <- list(list("d3xog_wt_nuc_478_rmdup", "d3xog_ko_nuc_256_rmdup"))
    rows <- 1
    columns <- 1
    orient <- 2
  }

  if (is.null(fname))  {x11("", 10, 5)
  } else {
    if (fname!="manual") {
      pdf(file=paste(profile2.path, "norm", "plots", fname, sep="/"), 6 * columns, 4.5 * rows)
    }  
  }                      
  par(mfrow=c(rows, columns), mar=c(2,4,1,1) + 0.1, oma=c(1, 1, 1, 1))
  if (orient==1) {
    data <- lapply(samples,
                 function(sample) lapply(sample, function(s)
                                         profileRead(paste(profile2.path, s, "profiles", sep="/"),
                                                           paste(annotation, data_type, sep="_"), group2)))
  } else if (orient==2) {
    data <- lapply(samples, function(sample)
                 lapply(sample, function(s) 
                        profileRead(paste(profile2.path, "norm",
                                          data_type, annotation, "profiles", sep="/"), s, group2)))
  }  
  ## If baseline is specified, normalize by mean start and end valuee
  #return(data)
  if (baseline) {
    data <- list(baselineNorm(data[[1]]))
  }  
#  return(data)
  
  if (!is.null(range)) {
      data <- list(lapply(data[[1]], function(x) lapply(x, function(y) as.matrix(y[range[1]:range[2],]))))
  }
 # return(data)
  if (is.null(y.vals)) {
    y.vals <- lapply(data, getRange)
    if (standard) {
      tmp_vals <- findMaxMin(y.vals)
      y.vals <- list(tmp_vals)
      for (i in 2:length(data)) y.vals <- c(y.vals, list(tmp_vals))
    }  
  } else {
    tmp_vals <- y.vals
    y.vals <- list(tmp_vals)
    for(i in 2:length(data)) y.vals <- c(y.vals, list(tmp_vals))
  }
  
  y_label_pos <- c(0.75, 0.25)
  lapply(c(1:length(samples)), function(x) {
    plotAnno(data[[x]], annotation, cols=cols, lab=lab, y.val=y.vals[[x]], wsize=wsize, stack=TRUE)
    #mtext("Normalized read count", at=y_label_pos[x], side=2, outer=TRUE, line=-1, cex=1)
    name <- lapply(samples[[x]], function(y) unlist(str_split(y[1], "_")))
    name <- unlist(lapply(name, function(y) y[-length(y)]))
    #legend(0, y=(y.vals[[x]][2] - (y.vals[[x]][2] / 20)), legend=legend, col=cols, bty="n", lty=1, horiz=TRUE)
  })
  par(las=1)
  mtext(annotation, side=3, outer=TRUE, cex=.8)
  if (!is.null(fname))
    if (fname!="manual") dev.off()
}

MP.plot2.horiz <- function(annotation, set="d3a", data_type = "raw", group2=NULL, cols=NULL, lab=c("",""), y.vals=NULL, wsize=25, standard=TRUE, fname=NULL) {
  orient <- 2
  if (set == "cells_pair") {
    samples <- list(list("icam_hmedip.bed", "icam_medip.bed"),
                    list("ngn_hmedip.bed", "ngn_medip.bed"),
                    list("omp_hmedip.bed", "omp_medip.bed"))
    rows <- 1
    columns <- 3
  } else if (set == "cells_rlm") {
   samples <- list(list("icam_hmc_rlm", "icam_mc_rlm"),                   
                   list("ngn_hmc_rlm", "ngn_mc_rlm"),
                   list("omp_hmc_rlm","omp_mc_rlm"))
    rows <- 1
    columns <- 3
  }  else if (set == "cells_norm") {
    samples <- list(list("icam_hmc", "icam_mc"),
                    list("ngn_hmc", "ngn_mc"),
                    list("omp_hmc", "omp_mc"))
    rows <- 1
    columns <- 3
  } else if (set=="cells_norm_hmc") {
    samples <- list("icam_hmc", "ngn_hmc", "omp_hmc")
    columns <- 3
    rows <- 1
    orient <- 2
  } else if (set=="cells_norm_mc") {
    samples <- list("icam_mc", "ngn_mc", "omp_mc")
    columns <- 3
    rows <- 1
    orient <- 2
    
  } else if (set == "d3a_pair") {
    samples <- list(list("moe_wt_hmc.bed"))
  } else if (set=="cells_rf") {
    samples <- list(list("omp_hmc_rf", "omp_mc_rf"))
    rows <- 1
    columns <- 1
    orient <- 1
  }
  
  if (is.null(fname))  {x11("", 10, 5)
  } else {
    pdf(file=paste(profile2.path, "plots", fname, sep="/"), 3 * columns, 3 * rows)
  }                      
  par(mfrow=c(rows, columns), mar=c(2,4,1,1) + 0.1, oma=c(1, 5, 1, 1))
  if (orient==1) {
    data <- lapply(samples, function(sample)
                 lapply(sample, function(s)
                        profileRead(paste(profile2.path, "norm", data_type, annotation, "profiles", sep="/"),
                                    s, group2)))
    
  } else if (orient==2) {
    data <- lapply(samples, function(sample)
                 lapply(sample, function(s)
                        profileRead(paste(profile2.path, "norm", data_type, annotation, "profiles", sep="/"),
                                    s, group2)))
  }  
            
  #return(data)

  data <- lapply(data, function(x) trimData(x, c(1, length(x[[1]][[1]]) - 1)))
  #return(data)
  if (is.null(y.vals)) {
    y.vals <- lapply(data, getRange)
    #y.vals <- lapply(data, function(x) lapply(x, getRange))
    print(y.vals)
    if (standard) {
      tmp_vals <- findMaxMin(y.vals)
      y.vals <- list(tmp_vals)
      for (i in 2:length(data)) y.vals <- c(y.vals, list(tmp_vals))
    }
    print(y.vals)
  } else {
    tmp_vals <- y.vals
    y.vals <- list(tmp_vals)
    for(i in 2:length(data)) y.vals <- c(y.vals, list(tmp_vals))
  }  
  lapply(c(1:length(samples)), function(x) {
    plotAnno(data[[x]], annotation, cols=cols, lab=lab, y.val=y.vals[[x]], wsize=wsize, stack=TRUE)
    #legend(0, y=(y.vals[[x]][2] - (y.vals[[x]][2] / 20)), legend=legend, col=cols, bty="n", lty=1, horiz=TRUE)
  })
  mtext("Normalized read count", at=.5, side=2, outer=TRUE, line=-1, cex=1)
  par(las=1)
  mtext(annotation, side=3, outer=TRUE, cex=1.6)
  if (!is.null(fname)) dev.off()
}



makeLegend <- function(labels, title=NULL, cols, fname=NULL) {
  if (is.null(fname)) {
        x11("",10,10)
      } else {
        pdf(paste(plot.path, paste(fname, ".pdf", sep=""), sep="/"), width=7, height=7)
  }
  plot(1, 1, type="n", xlim=c(0,10), ylim=c(0,10))
  legend(3, 5, legend=labels, title=title, col=cols, lwd=2, bty="n", lty=1, horiz=TRUE, cex=1, adj=c(1.2,-2), x.intersp=.6)
  if (!is.null(fname)) dev.off()
    
}

prepMP <- function(vals) {
  vals <- vals[is.finite(vals$ams_A),]
  vals <- vals[!is.na(vals$ams_A),]
  return(vals)
}

subsetByROI <- function(data, roi, select=4, norm_factor=1) {
  raw <- vector(length=nrow(roi), mode="numeric")
  norm <- vector(length=nrow(roi), mode="numeric")
  #thresh <- roi[1,3] - roi[1,2]
  #x <- 1
  i <- 1
  name <- ""
  #while (i <= nrow(roi)) {
  for (j in 1:nrow(data)) {
    if (data[j,2] >= roi[i,2]) {
      dist_j1 <- abs(data[j-1,2] - roi[i,2])
      dist_j2 <- abs(data[j,2] - roi[i,2])
      ind <- 0
      if (dist_j1 < dist_j2) {
        ind <- j-1 
      } else {
        ind <- j
      }
      raw[i] <- data[ind, 3]
      norm[i] <- data[ind, 4]
      if (i < nrow(roi)) {
        i <- i + 1
      } else {
        break
      }  
        #x <- j + 1
        #break
    }
  }
  raw <- raw / norm_factor
  out <- cbind(roi, raw, norm)
  return(out)
}    
 
subsetByROI.par <- function(sample, anno.path, roi_name, select=4) {
  roi <- read.delim(paste(anno.path, roi_name, sep="/"), header=FALSE)
  roi <- roi[order(roi[,1], roi[,2]),]
  random_ind <- grep("random", roi[,1])
  if (length(random_ind) > 0) roi <- roi[-random_ind,]
  chr_names <- unique(roi[,1])
  roi_split <- split(roi, roi[,1])
  
  norm_factor <- scan(file=paste(medip_split.path, sample, "total_reads", sep="/")) / 1e6
  registerDoMC(cores=6)
  out <- foreach (chr=chr_names, .combine="rbind") %dopar% {
    cat("=")
    data.path <- paste(medip_split.path, sample, chr, sep="/")
    data <- read.delim(data.path, header=FALSE)
    #print(as.character(chr))
    return(subsetByROI(data, roi_split[[chr]], select, norm_factor))
  }
  cat("\n")
  #return(out)
  colnames(out) <- c("chr", "start", "end", "name", "group", "strand", "raw", "norm")
  out.path <- paste(profile2.path, sample, sep="/")
  if (!file.exists(out.path)) dir.create(out.path)
  cat("Writing...\n")
  write.table(out, file=paste(out.path, roi_name, sep="/"), quote=FALSE, sep="\t", row.names=FALSE)
}


saveRoiByChr <- function(roi_name) {
  roi <- read.delim(roi_name, header=F)
  rand <- grep("random", roi[,1])
  if (length(rand) > 0) roi <- roi[-rand,]
  roi.split <- split(roi, roi[,1])
  #return(roi.split)
  #rand <- grep("random", names(roi.split))
  #no_rand <- c(1:length(roi.split))[-rand]
  #print(rand)
  #print(no_rand)
  #if (length(rand) > 0) roi.split <- roi.split[[no_rand]]
  roi.split <- lapply(roi.split, function(x) x[order(x[,1], x[,2]),])
  chrs <- names(roi.split)
  roi_split_name <- paste(roi_name, "_chr", sep="")
  if (!file.exists(roi_split_name)) {
    dir.create(roi_split_name)
  } else {
    return
  }
  sapply(chrs, function(chr) {
    write.table(roi.split[[chr]], file=paste(roi_split_name, chr, sep="/"),
                quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
  })
}

saveRoiByChr.all <- function(anno_set="std") {
  if (anno_set=="std") {
    anno.path <- anno.path
  } else if (anno_set == "hires") {
    anno.path <- anno_hires.path
  } else if (anno_set == "features") {
    print("here")
    anno.path <- feature.path
  }
  files <- list.files(anno.path)
  ind <- grep("chr", files)
  if (length(ind) > 0) files <- files[-ind]
  a <- foreach(file=files) %dopar% {
    if (!file.exists(paste(anno.path, paste(file, "_chr", sep=""), sep="/"))) {
      cat(file)
      cat("\n")
      result <- try(saveRoiByChr(paste(anno.path, file, sep="/")))
      if (class(result) == "try-error") return
    }
}
}

