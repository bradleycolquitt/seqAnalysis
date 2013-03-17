suppressPackageStartupMessages(library(itertools))
suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(boot))
library(RColorBrewer)
suppressPackageStartupMessages(library(multicore))
suppressPackageStartupMessages(library(doMC))
library(MEDIPS)
suppressPackageStartupMessages(library(colorspace))
suppressPackageStartupMessages(library(cluster))
suppressPackageStartupMessages(library(gplots))
library(stringr)

source("~/src/seqAnalysis/R/paths.R")
source("~/src/seqAnalysis/R/boot.R")
source("~/src/seqAnalysis/R/plotUtil.R")
source("~/src/seqAnalysis/R/samples.R")
source("~/src/MEDIPS/R/MEDIPS_mod.methylProfiling.R")

registerDoMC(cores=2)

col4 <- brewer.pal(4, "Spectral")
col3 <- brewer.pal(3, "Dark2")
col2 <- rev(brewer.pal(3, "Set1")[1:2])

# Make Profile ------------------------------------------------------------


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
  print(out.name)
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
        print(paste(out.path, paste(out.name, name, fun, "bootCI", x, sep="_")))
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
    return(tapply(data[,param[[1]]], grouping, param[[3]], na.rm=TRUE))
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
makeProfile2 <- function(anno, samples, group2=NULL, data_type="unnorm/mean", fun="mean", rm.outliers=0, 
                         sample_num=0, exclude_chr=NULL, write=TRUE) {
  sample_path <- paste(profile2.path, "norm", data_type, anno, sep="/")
  a <- sapply(samples, function(sample) {
      print(sample)
      data <- read.delim(paste(sample_path, sample, sep="/"))
      colnames(data) <- c("chr", "start", "end", "name", "group", "strand", "norm")
      ind <- "norm"
      profile <- NULL
      CI <- NULL
      
      if (!is.null(exclude_chr)) {
        data <- data[!(data$chr %in% exclude_chr),]
      }
      
      if (sample_num > 0) {
        print(sample_num)
        sample_names <- sample(unique(data$name), sample_num)
        data <- data[data$name %in% sample_names,]
        gc()
      }

      if (rm.outliers > 0) {
        thresh <- quantile(data$norm, probs=c(rm.outliers, 1-rm.outliers), na.rm=TRUE)
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
        if (!is.null(exclude_chr)) {
          sample <- paste(sample, paste("exclude", exclude_chr, sep="-"), sep="_")
          
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
makeProfile2.allSamp <- function(anno, group2=NULL, data_type="unnorm/mean", fun="mean", rm.outliers=0, sample_num=0, exclude_chr=NULL, write=T) {  
  sample_path <- paste(profile2.path, "norm", data_type, anno, sep="/")
  print(sample_path)
  samples <- list.files(sample_path)
  ind <- grep("profiles", samples)
  if (length(ind) > 0) samples <- samples[-ind]
  out_path = paste(sample_path, "profiles", sep="/")
  data <- foreach(sample=samples) %do% {
    out_name <- sample
    if (!is.null(exclude_chr)) {
      out_name <- paste(out_name, paste("exclude", exclude_chr, sep="-"), sep="_")
    }
    if (!is.null(group2)) {
      out_name <- paste(out_name, group2, sep="_")
    }
    if (rm.outliers > 0) {
      out_name <- paste(out_name, paste("trim", rm.outliers, sep=""), sep="_")
    }
    if (file.exists(paste(out_path, paste(out_name, fun, sep="_"), sep="/"))) {
      print("Skipping")
      next 
    }
    print(out_name)
    return(makeProfile2(anno, sample, data_type=data_type, fun=fun, group2=group2,
                           rm.outliers=rm.outliers, sample_num=sample_num, exclude_chr=exclude_chr, write=write))
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


# Plot Profiles -----------------------------------------------------------


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
       ann=FALSE, axes=FALSE, ...)

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
  
  # Draw and label axes
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

## Primary interface for drawing profile of a single sample 
plot2 <- function(annotation, sample, orient=2, data_type="unnorm/mean", fun="mean", group2=NULL, group2_col=NULL,
                     cols=1, lab=c("",""), y.vals=NULL, wsize=25, range=NULL, type="range", fname=NULL, ...) {

  ## Orient 1 data structure is no longer used (2/14/12)
  if (orient==1) {
    data <- profileRead(paste(profile2.path, sample, "profiles", sep="/"), fun,
                        paste(annotation, group2, data_type, sep="_"))
  } else if (orient==2) {
    data <- profileRead(paste(profile2.path, "norm", data_type, annotation, "profiles", sep="/"), fun, sample, group2)
  }
 # return(data)
  ## If group2_col is specified, select given column of data
  if (!is.null(group2_col)) {
    if (length(group2_col) == 1) {
      data <- lapply(data, function(x) matrix(x[,group2_col], ncol=1))
    } else {
      data <- lapply(data, function(x) x[,group2_col])
    }
  }
  #return(data)
  ## If range specificed, trim profile by given values
  if (!is.null(range)) {
      data <- lapply(data, function(x) matrix(x[range[1]:range[2],], ncol=1))
  }
  #return(data)
  ## Set up graphics device
  if (is.null(fname))  {
    x11("", 6, 4)
  } else if (fname=="manual") {
    # do nothing, device has already been set
  } else {
    if (fname=="auto") {
      dt <- paste(unlist(str_split(data_type, "/")), collapse="_")
      if (!is.null(group2)) {
        fname <- paste(annotation, sample, dt, group2, fun, sep="_")
      } else {
        fname <- paste(annotation, sample, dt, fun, sep="_")        
      }
    }  
    fname <- paste(fname, ".pdf", sep="")
    print(paste("Saving to ", fname, sep=""))
    pdf(file=paste(profile2.path, "norm", "plots", fname, sep="/"), 6, 4.5)
  }
  
  ## If y.vals not provided, determine y axis range from the data
  if (is.null(y.vals)) y.vals <- getRange(list(data), buffer=0)
  
  ## Send data to plotAnno
  plotAnno(list(data), annotation, cols=cols, lab=lab,
                y.val=y.vals, wsize=wsize, ...) 
  
  mtext(paste(sample, " -> ", group2, sep=""), side=3, outer=T, cex=1.6)
  if (!is.null(fname)) {
    if (fname != "manual") dev.off()
  }
}

## Plot several profiles in one graphic device
plot2.several <- function(annotation, set="d3a", data_type="unnorm/mean", group2=NULL, cols=NULL, lab=c("",""), 
                          y.vals=NULL, wsize=25, standard=FALSE, range=NULL, baseline=FALSE, group2_col=NULL,
                          fun="mean", legend=FALSE, fname=NULL) {
  samples <- NULL
  orient <- 1
  rows <- 1
  columns <- 1

  plot_param <- select_set(set)
  samples <- plot_param[[1]]
  rows <- plot_param[[2]]
  columns <- plot_param[[3]]
  orient <- plot_param[[4]]

  # setup graphic device
  if (is.null(fname))  {
    x11("", 5, 6)
  } else if (fname=="manual") {
    # do nothing, device has already been set
    
  } else {
    if (fname=="auto") {
      dt <- paste(unlist(str_split(data_type, "/")), collapse="_")
      if (!is.null(group2)) {
        fname <- paste(annotation, set, dt, group2, fun, paste("y", paste(y.vals, collapse="_"), sep=""), sep="_")
      } else {
        fname <- paste(annotation, set, dt, fun, paste("y", paste(y.vals, collapse="_"), sep=""), sep="_")        
      }
    }  
    fname <- paste(fname, ".pdf", sep="")
    print(paste("Saving to ", fname, sep=""))
    pdf(file=paste(profile2.path, "norm", "plots", fname, sep="/"), 6, 9)
  }
  par(mfrow=c(rows, columns))
      #mar=c(2,4,1,1) + 0.1, oma=c(1, 1, 1, 1))
  
  # Read in data
  if (orient==1) {
    data <- lapply(samples,
                 function(sample) lapply(sample, function(s)
                                         profileRead(paste(profile2.path, s, "profiles", sep="/"),
                                                           paste(annotation, data_type, sep="_"), group2)))
  } else if (orient==2) {
    data <- lapply(samples, function(sample)
                 lapply(sample, function(s) 
                        profileRead(paste(profile2.path, "norm",
                                          data_type, annotation, "profiles", sep="/"), fun, s, group2)))
  }
  #return(data)
  ## If group2_col is specified, select given column of data
  if (!is.null(group2_col)) {
    if (length(group2_col) == 1) {
      data <- lapply(data, function(x) lapply(x, function(y) lapply(y, function(z) matrix(z[,group2_col], ncol=1))))
    } else {
      data <- lapply(data, function(x) lapply(x, function(y) lapply(y, function(z) z[,group2_col])))
    }
  }
  #if (!is.null(group2_col)) {
  #  if (length(group2_col) == 1) {
  #    data <- lapply(data, function(x) matrix(x[,group2_col], ncol=1))
  #  } else {
  #    data <- lapply(data, function(x) x[,group2_col])
  #  }
  #}
  ## If baseline is specified, normalize by mean start and end values  
  if (baseline) {
    data <- list(baselineNorm(data[[1]]))
  }  

  ## Trim data to x-axis range specificed with 'range' parameter
  if (!is.null(range)) {
      data <- list(lapply(data[[1]], function(x) lapply(x, function(y) as.matrix(y[range[1]:range[2],]))))
  }
  
  ## If y.vals unspecficied, calculate reasonable range from the data
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
  
  ## Send data for each sample to the plotter
  lapply(c(1:length(samples)), function(x) {
    plotAnno(data[[x]], annotation, cols=cols, lab=lab, y.val=y.vals[[x]], wsize=wsize, stack=TRUE)
    name <- lapply(samples[[x]], function(y) unlist(str_split(y[1], "_")))
    name <- unlist(lapply(name, function(y) y[-length(y)]))
    if (legend) {
      legend(2, y.vals[[x]][2], legend=samples[[x]], col=cols, lty=1)
    }
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
  if (is.null(fname))  {
    x11("", 6, 4)
  } else if (fname=="manual") {
    # do nothing, device has already been set
  } else {
    if (fname=="auto") {
      dt <- paste(unlist(str_split(data_type, "/")), collapse="_")
      if (!is.null(group2)) {
        fname <- paste(annotation, set, dt, group2, fun, sep="_")
      } else {
        fname <- paste(annotation, set, dt, fun, sep="_")        
      }
    }  
    fname <- paste(fname, ".pdf", sep="")
    print(paste("Saving to ", fname, sep=""))
    pdf(file=paste(profile2.path, "norm", "plots", fname, sep="/"), 6, 4.5)
  }
  
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



saveRoiByChr <- function(roi_name) {
  roi <- read.delim(roi_name, header=F)
  rand <- grep("random", roi[,1])
  if (length(rand) > 0) roi <- roi[-rand,]
  roi.split <- split(roi, roi[,1])
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

