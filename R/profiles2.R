
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

source("~/src/R/paths.R")
source("~/src/R/boot.R")
source("~/src/R/plotUtil.R")
source("~/src/MEDIPS/R/MEDIPS_mod.methylProfiling.R")

registerDoMC(cores=10)

col4 <- brewer.pal(4, "Spectral")
col3 <- brewer.pal(3, "Dark2")
col2 <- rev(brewer.pal(3, "Set1")[1:2])

writeModule <- function(out.path, sample, group2=NULL, profile, CI) {
  if (!file.exists(out.path)) {
    dir.create(out.path)
  }
  out.name <- NULL

  if (is.null(group2)) {
    out.name <- sample
  } else {
    out.name <- paste(sample, group2, sep="_")
  }
  if (length(profile) > 1) {
    lapply(names(profile), function(name) {
      write.table(as.matrix(profile[[name]]),
                  file=paste(out.path, paste(out.name, name, "mean", sep="_"), sep="/"),
                  quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
      lapply(c(1:2), function(x) {
        unwrapped <- NULL
        if (is.null(group2)) {unwrapped <- .unwrap(CI[[name]], x)
        } else {unwrapped <- apply(CI[[name]], 2, function(y) .unwrap(y, x))}
        write.table(unwrapped,
                  file=paste(out.path, paste(out.name, name, "mean_bootCI", x, sep="_"),
                    sep="/"),
                    quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
      })
    })
    
  } else {
    write.table(profile,
              file=paste(out.path, paste(out.name, "mean", sep="_"), sep="/"),
              quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)     
    
    lapply(c(1:2), function(x) {
      unwrapped <- NULL
      if (is.null(group2)) {unwrapped <- .unwrap(CI[[1]], x)
      } else {unwrapped <- apply(CI[[1]], 2, function(y) .unwrap(y, x))}
      write.table(unwrapped,
                  file=paste(out.path, paste(out.name, "mean_bootCI", x, sep="_"),
                    sep="/"),
                    quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
    })
  }  
}


profileCompute <- function(data, param) {
  grouping <- lapply(param[[2]], function(x) if(!is.null(x)) data[,x])
  return(tapply(data[,param[[1]]], grouping, param[[3]]))
}

MP.makeProfile2 <- function(anno, sample, group2=NULL, data_type="rpm_avg", write=TRUE) {
  data <- read.delim(paste(anno, sample, sep="/"))
  colnames(data) <- c("chr", "start", "end", "name", "group", "strand", "norm")
  #ind <- c("raw", "norm")
  ind <- "norm"
  profile <- NULL
  CI <- NULL
#  return(data)
  if (!is.null(group2)) {
    group2_vals <- read.delim(paste(group2.path, group2, sep="/"), header=FALSE)
    data$group2 <- group2_vals[match(data$name, group2_vals[,1]),2]
    group2_name <- "group2"
    profile <- lapply(ind, function(x) profileCompute(data, list(x, list("group", group2_name), mean)))
    names(profile) <- ind
    CI <- lapply(ind, function(x) profileCompute(data, list(x, list("group", group2_name), bootCI)))
    names(CI) <- ind
  } else {
    profile <- lapply(ind, function(x) profileCompute(data, list(x, list("group"), mean)))
    names(profile) <- ind
    CI <- lapply(ind, function(x) profileCompute(data, list(x, list("group"), bootCI)))
    names(CI) <- ind
  }
  #return(CI) 
  if(!write) {
    return(profile)
  } else {
    out.path <- paste(anno, "profiles", sep="/")
    if (!is.null(group2)) {
      group2_out <- str_split(group2, "/")
      group2_out <- group2_out[length(group2_out)]
    } else {
      group2_out <- NULL
    }   
    writeModule(out.path=out.path, sample=sample, group2=group2_out, profile=profile, CI=CI)   
  }  
}

MP.makeProfile2.allSamp <- function(anno, set="d3a", group2=NULL, data_type="rpm_avg", write=T) {
  if (set=="d3a") {samples <- samples.d3a
  } else if (set=="cells") {samples <- samples.cells
  } else if (set=="cells_rlm") {samples <- samples.cells.rlm}
  #samples <- lapply(samples, function(sample) paste(profile2.path, sample, sep="/"))
  sample_path <- paste(profile2.path, "norm", data_type, anno, sep="/")
  print(sample_path)
  samples <- list.files(sample_path)
  print(samples)
  ind <- grep("profiles", samples)
  if (length(ind) > 0) samples <- samples[-ind]
  out_path = paste(sample_path, "profiles", sep="/")
  #if (norm) samples <- paste("norm", samples, sep="/")
  data <- foreach(sample=samples) %dopar% {
    if (file.exists(paste(out_path, sample, sep="/"))) {
      return
    }
    print(sample)
    return(MP.makeProfile2(sample_path, sample, group2=group2, write=write))
  }
  names(data) <- samples
  return(data)
}

MP.makeProfile2.allAnno <- function(sample, set="d3a", group2=NULL, write=TRUE) {
  files <- list.files(paste(profile2.path, "norm", sep="/"))
  registerDoMC(cores=4)
  i <- 0
  for(file in files) {
    i <- i + 1
    cat("--")
    cat(file)
    cat(paste(" ", i, " of ", length(files), sep=""))
    cat("\n")
    tryCatch(MP.makeProfile2.allSamp(file, group2=group2, rpm=TRUE, write=write), error = function(e) {
      print(paste("Skipping", file, sep=" "))
      print(e)
      return
    })
#    }
  } 
}

MP.makeProfile2.all <- function(set="d3a", write=TRUE) {
  if (set=="d3a") samples <- samples.d3a
  if (set=="cells") samples <- samples.cells
  if (set=="cells_rlm") samples <- samples.cells.rlm
  if (set=="d3a_rlm") samples <- samples.d3a.rlm
  for (sample in samples) {
    cat(sample)
    cat("\n")
    sample.path <- paste(profile2.path, sample, sep="/") 
    MP.makeProfile2.allAnno(sample.path, write=write)
  }
}
MP.makeProfile.AllAnno <- function(sample, group2=NULL, select=NULL, fname) {
  file.path <- paste("~/storage/analysis/mprofiles", sample, sep="/")
  profile.path <- paste(file.path, "profiles", sep="/")
  if (!file.exists(profile.path)) dir.create(profile.path)
  extend <- NULL
  if (!is.null(group2)) extend <- group2
  if (!is.null(select)) extend <- select
  files.tbp <- checkExisting(file.path, profile.path, extend=extend)  
  files.tbp <- files.tbp[-grep("profiles", files.tbp)]
  for(file in files.tbp) {
    cat(file)
    cat("\n")
    result <- try(MP.makeProfile(sample=sample, mp=file, group2=group2, select=select, fname=fname))
    if (class(result) == "try-error") next
  }
}

MP.makeProfile.AllSamples <- function(mp, group2=NULL, select=NULL, cluster=NULL, fname) {
  x <- foreach(sample=samples) %dopar% {
    file.path <- paste("~/storage/analysis/mprofiles", sample, sep="/")
    profile.path <- paste(file.path, "profiles", sep="/")
    if (!file.exists(profile.path)) dir.create(profile.path)
    extend <- NULL
    if (!is.null(group2)) extend <- group2
    if (!is.null(select)) extend <- select
    cat(sample)
    cat("\n")
    MP.makeProfile(sample=sample, mp=mp, group2=group2, select=select, cluster=cluster, fname=fname)
  }
}

MP.makeProfile.Genes <- function(sample, group2=NULL) {
  file.path <- paste("~/analysis/mprofiles", sample, sep="/")
  profile.path <- paste(file.path, "profiles", sep="/")
  if (!file.exists(profile.path)) dir.create(profile.path)
  files.tbp <- checkExisting(file.path, profile.path, select=gene.profiles, extend=group2)
  for(file in files.tbp) {
    cat(file)
    cat("\n")
    MP.makeProfile(sample, file, group2)
  }
}

MP.positionMatrix <- function(data, data_type="raw") {
  groups <- unique(data$group)
  names <- unique(data$name)
  #out <- matrix(NA, nrow=length(names), ncol=length(groups), dimnames=list(names, groups))
  data_groups = split(data, data$group)
  
  out <- foreach(ind=icount(length(data_groups)), .combine="cbind") %dopar% {
    data_group <- data_groups[[ind]]
    m <- match(names, data_group$name)
    return(data_group[m, data_type])
    #out[,ind] <- data_group[m, data_type]
  }
  rownames(out) <- names
  #for(i in 1:length(groups)) {
  #  tmp <- subset(mp, mp$group==groups[i])
  #  m <- match(rownames(out), tmp$name)
  #  out[,i] <- tmp[m, data_type]
  #}
  return(out)
}

MP.positionMatrix.all <- function(anno, set="d3a", data_type="raw") {
  if (set=="d3a") {
    samples <- samples.d3a
  } else if (set=="cells") {
    samples <- samples.cells
  }
  for (sample in samples) {
    print(sample)
    data <- read.delim(paste(profile2.path, sample, anno, sep="/"), header=FALSE, colClasses=profile.classes)
    colnames(data) <- c("chr", "start", "end", "name", "group", "strand", "raw", "norm")
    pos_matrix <- MP.positionMatrix(data, data_type=data_type)
    out_path <- paste(profile2.path, sample, "images", sep="/")
    if (!file.exists(out_path)) dir.create(out_path)
    write.table(pos_matrix, file=paste(profile2.path, sample, "images", anno, sep="/"),
                quote=FALSE, sep="\t", col.names=FALSE)
  }
  
}


MP.makeImage <- function(sample, anno, image=TRUE) {
  image.path <- paste(profile2.path, sample, "images", anno, sep="/")
  print(image.path)
  data <- read.delim(image.path, header=FALSE, row.names=1)
  #mp.data <- prepMP(mp.data)
  #mp.data <- thresh
  #out <- MP.positionMatrix(mp.data)
  if (image) {MP.image(data)}
  return(data)
}

MP.makeImages.all <- function(set="d3a", anno) {
  if (set=="d3a") {
    samples <- samples.d3a
  }
  out <- lapply(samples, function(sample) {
    MP.makeImage(sample, anno, image=FALSE)
  }) 

  return(out)
  
}

MP.processImage <- function(vals) {
  v <- apply(vals, 2, pseudoCountNorm)
  v <- na.omit(v)
  v <- log(v, 2)
  return(v)
}
MP.image <- function(vals) {
  X11()
  image(t(vals),
        #col=rev(heat_hcl(50, c=c(80,30), l=c(30,90), power=c(1/5, 1.3))))
        col=greenred(50))
  
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

MP.plotProfiles <- function(profile=NULL, ci_1=NULL, ci_2=NULL, cols=NULL,
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
  
  if (!is.null(ci_1)) {
    #if (!is.null(dim(profile)))
      polygon(c(c(1:length(profile)), c(length(profile):1)),
              c(ci_1, rev(ci_2)),
              col=rgb(col.rgb[1,1], col.rgb[2,1], col.rgb[3,1],
                alpha=col.rgb[4,1], maxColorValue=255),
              border=F)
  }
  lines(profile, col=cols)
}

MP.plotAnno <- function(data, annotation, cols=NULL, lab=NULL,
                        y.val=NULL, combine=FALSE,
                        stack=FALSE, ...) {
  lab.data <- NULL
  if (is.null(dim(data[[1]][[1]]))) {
    x.lim <- length(data[[1]][[1]])
  } else {
    x.lim <- nrow(data[[1]][[1]])
  }
 
  ## Send data to line drawer
  #if (stack || length(data) > 1) {
    plot(1, 1, type="n",
       xlim=c(1, x.lim),
       ylim=c(round(y.val[1],2), round(y.val[2], 2)),
       xlab="",
       ylab="",
       ann=FALSE, axes=FALSE)
    
    for(i in 1:length(data)) {
      data.curr <- data[[i]]
      data.val <- data.curr[[1]]
      data.ci1 <- data.curr[[2]]
      data.ci2 <- data.curr[[3]]
      if (ncol(data.val) > 1) {
        for (i in 1:ncol(data.val)) {
          MP.plotProfiles(profile=data.val[,i], ci_1=data.ci1[,i], ci_2=data.ci2[,i],
                          smooth=FALSE, cols=cols[i])
          
        }
        lab.data <- computeAxis(data.val[,1], lab)
      } else {

        MP.plotProfiles(profile=data.val[,1], ci_1=data.ci1[,1], ci_2=data.ci2[,1],
                          smooth=FALSE, cols=cols[i])
        lab.data <- computeAxis(data.val[,1], lab)
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
  abline(v=lab.data$pos[3], lty=2, col="grey")
}




MP.plot2 <- function(annotation, sample, data_type = "rpm_avg", group2=NULL, cols=NULL, lab=NULL, y.vals=NULL, type="range", fname=NULL) {
  ## Get sample names
  if (!is.null(group2)) {
    anno <- paste(annotation, group2, data_type, sep="_")
  } else {
    anno <-  paste(annotation, data_type, sep="_")
  }
 
  data <- profileRead(paste(profile2.path, "norm", data_type, annotation, "profiles", sep="/"), sample, group2)
  #data <- MP.getData(profile2.path, sample, annotation=anno, group2=group2)[[1]]
  #return(data)
  #data <- splitReform(data)
  #if (!is.null(dim(data[[1]]))) {
  #  #print("here")
  #  data <- lapply(data, function(x) apply(x, 2, function(y) y[1:length(y) - 1]))
  #} else {
  #  data <- lapply(data, function(x) trimData(x, c(0, length(x[[1]][[1]]) - 1)))
                                        #return(data)
  #}
  if (is.null(fname))  {x11("", 6, 4)
  } else {
    pdf(file=paste(profile2.path, "norm", "plots", fname, sep="/"), 6, 4.5)
  }
  .profile <- function(data, y.vasl) {
    MP.plotAnno(list(data), annotation, cols=cols, lab=lab,
                y.val=y.vals) 
  }
  if (is.null(y.vals)) y.vals <- getRange(list(data), buffer=0)
  #print(y.val)
  .profile(data, y.vals)
  
  #mtext("AMS", at=.775, side=2, outer=T, line=-1, cex=1.6)
  #mtext("AMS", at=.275, side=2, outer=T, line=-1, cex=1.6)
  mtext("Normalized read count", at=.5, side=2, outer=T, line=-2, cex=1)
  par(las=1)
  #mtext("5hmC", at=.775, side=2, outer=T, line=1, cex=1.6)
  #mtext("5mC", at=.275, side=2, outer=T, line=1, cex=1.6)
  #mtext("5hmC", at=.275, side=3, outer=T, line=0, cex=1.6)
  #mtext("5mC", at=.775, side=3, outer=T, line=0, cex=1.6)
  #mtext(paste(sample, " -> ", group2, sep=""), side=3, outer=T, cex=1.6)
  if (!is.null(fname)) {
    dev.off()
  }
}

MP.plot2.several <- function(annotation, set="d3a", data_type="rpm", group2=NULL, cols=NULL, lab=NULL, y.vals=NULL, standard=FALSE,
                             fname=NULL) {
  samples <- NULL
  orient <- 1
  if (set=="d3a") {
    samples <- list(list("moe_wt_hmc.bed", "moe_d3a_hmc.bed"), list("moe_wt_mc.bed", "moe_d3a_mc.bed"))
    legend <- c("Dnmt3a +/+", "Dnmt3a -/-")
    rows <- 2
    columns <- 1
  } else if (set=="d3a_rlm") {
    samples <- list(list("moe_wt_hmc_rlm", "moe_d3a_hmc_rlm"), list("moe_wt_mc_rlm", "moe_d3a_mc_rlm"))
    rows <- 2
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
    if (!is.null(group2)) {
      columns <- 3
    }
    rows <- 2
    
    orient <- 2
  }

  #print(samples)
  #columns <- length(samples[[1]])
  if (is.null(fname))  {x11("", 10, 5)
  } else {
    pdf(file=paste(profile2.path, "plots", fname, sep="/"), 6 * columns, 3 * rows)
  }                      
  
  par(mfrow=c(rows, columns), mar=c(2,4,1,1) + 0.1, oma=c(1, 5, 1, 1))
  if (orient==1) {
    data <- lapply(samples,
                 function(sample) lapply(sample, function(s)
                                         profileRead(paste(profile2.path, s, "profiles", sep="/"),
                                                           paste(annotation, data_type, sep="_"), group2)))
  } else if (orient==2) {
    data <- lapply(samples, function(sample)
                 lapply(sample, function(s)
                        profileRead(paste(profile2.path, "norm", data_type, annotation, "profiles", sep="/"), s, group2)))
  }  

  #return(data)

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
    MP.plotAnno(data[[x]], annotation, cols=cols, lab=lab, y.val=y.vals[[x]], stack=TRUE)
    mtext("Normalized read count", at=y_label_pos[x], side=2, outer=TRUE, line=-1, cex=1)
    name <- lapply(samples[[x]], function(y) unlist(str_split(y[1], "_")))
    name <- unlist(lapply(name, function(y) y[-length(y)]))
    #legend(0, y=(y.vals[[x]][2] - (y.vals[[x]][2] / 20)), legend=legend, col=cols, bty="n", lty=1, horiz=TRUE)
  })
  par(las=1)
  mtext("5hmC", at=.775, side=2, outer=T, line=1, cex=1.6)
  mtext("5mC", at=.275, side=2, outer=T, line=1, cex=1.6)
  mtext(annotation, side=3, outer=TRUE, cex=1.6)
  if (!is.null(fname)) dev.off()
}

MP.plot2.horiz <- function(annotation, set="d3a", data_type = "raw", group2=NULL, cols=NULL, lab=c("",""), y.vals=NULL, standard=TRUE, fname=NULL) {
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
  }
  else if (set == "d3a_pair") {
    samples <- list(list("moe_wt_hmc.bed"))
  }
  
  if (is.null(fname))  {x11("", 10, 5)
  } else {
    pdf(file=paste(profile2.path, "plots", fname, sep="/"), 3 * columns, 3 * rows)
  }                      
  par(mfrow=c(rows, columns), mar=c(2,4,1,1) + 0.1, oma=c(1, 5, 1, 1))
  data <- lapply(samples, function(sample)
                 lapply(sample, function(s)
                        profileRead(paste(profile2.path, "norm", annotation, "profiles", sep="/"), s, group2)))
            
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
    MP.plotAnno(data[[x]], annotation, cols=cols, lab=lab, y.val=y.vals[[x]], stack=TRUE)
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
    cat(file)
    cat("\n")
    result <- try(saveRoiByChr(paste(anno.path, file, sep="/")))
    if (class(result) == "try-error") next
  }
}


