
library(itertools)
library(foreach)
library(boot)
library(RColorBrewer)

col4 <- brewer.pal(4, "Spectral")
col3 <- brewer.pal(3, "Dark2")

AnnoSetSimp.makeProfile <- function(data=NULL, reads=NULL, data_type="rms",
                                    fname=NULL, write=TRUE, boot=FALSE, FUN=mean) {
  group <- data@set_group

  data_name <- unlist(strsplit(data@name, "/"))
  data_name <- data_name[length(data_name)]
  out_path <- paste("~/analysis/profiles", data_name, sep="/")
  #if (!file.exists(out_path)) {dir.create(out_path)}
  reads <- .extract(data@read_data, reads)
  norm_reads <- getDataByType(reads, data_type)
  
  if (boot) FUN <- AnnoSetSimp.bootCImean
  profile <- tapply(norm_reads, group, FUN)
  if (write) {
    if (!file.exists(out_path)) {
      dir.create(out_path)
    }
    if (boot) {
      lapply(c(1:2), function(x) {
        unwrapped <- .unwrap(profile, x) 
        
        write.table(unwrapped, file=paste(out_path, paste(fname, "bootCI", x, sep="_"), sep="/"),
                    quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
      })
    } else {
      write.table(as.matrix(profile), file=paste(out_path, fname, sep="/"), sep="\t", quote=FALSE,
                row.names=FALSE, col.names=FALSE, )
    }
  }
  if (!write) {return(profile)}
}

AnnoSetSimp.makeProfile.All <- function(data=NULL, data_type=norm, write=TRUE,
                                        FUN=mean, boot=FALSE, funcName="mean") {
  read_names <- .getNames(data@read_data)
  x <- foreach(curr_read_data=read_names) %dopar% {
    data_name <- unlist(strsplit(curr_read_data, ".bed"))[1]
    fname <- paste(data_name, funcName, sep="_")
    AnnoSetSimp.makeProfile(data, reads=curr_read_data, data_type=data_type,
                                 write=write, fname=fname, FUN=FUN, boot=boot)
  }
}


AnnoSetSimp.makeProfileGroup <- function(data=NULL, reads=NULL, data_type=norm, group2=NULL,
                                         fname=NULL, write=TRUE, boot=FALSE, FUN=mean) {
  group <- data@set_group
  group2_path <- "~/analysis/rna/quantiles"
  data_name <- unlist(strsplit(data@name, "/"))
  data_name <- data_name[length(data_name)]
  out_path <- paste("~/analysis/profiles", data_name, sep="/")
  group2 <- read.delim(paste(group2_path, group2, sep="/"), header=FALSE)
  group2 <- group2[group2[,2]>0,]
  names <- data@set_names
  reads <- .extract(data@read_data, reads)
  reads_vals <- getDataByType(reads, data_type)
  group2_vals <- group2[match(names, group2[,1]), 2]

  if (boot) FUN <- AnnoSetSimp.bootCImean
  profile <- tapply(reads_vals, list(group, group2_vals), FUN)
  if (write) {
    if (!file.exists(out_path)) {
      dir.create(out_path)
    }
    if (boot) {
      lapply(c(1:2), function(x) {
        unwrapped <- apply(profile, 2, function(y) {
          out <- .unwrap(y, x)
          return(out)
        })
        write.table(unwrapped, file=paste(out_path, paste(fname, "bootCI", x, sep="_"), sep="/"),
                    quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
      })
    } else {
      write.table(as.matrix(profile), file=paste(out_path, fname, sep="/"), sep="\t", quote=FALSE,
                row.names=FALSE, col.names=FALSE, )
    }
  }  
  return()
}

AnnoSetSimp.makeProfileGroup.All <- function(data, data_type="norm", group2=NULL, write=TRUE,
                                             FUN=mean, boot=FALSE, funcName="mean") {
  read_names <- .getNames(data@read_data)
  x <- foreach(curr_read_data=read_names) %dopar% {
    data_name <- unlist(strsplit(curr_read_data, ".bed"))[1]
    fname <- paste(data_name, group2, data_type, funcName, sep="_")
    AnnoSetSimp.makeProfileGroup(data, reads=curr_read_data, data_type=data_type,
                                 group2=group2, write=write, fname=fname, FUN=FUN, boot=boot)
  }
}


plotProfiles.TSS <- function(profile_matrix, ...) {
  x11("", 7, 7)
  plot(1, 1, type="n",
       xlim=c(1,nrow(profile_matrix)),
       ylim=c((min(profile_matrix)-10), (max(profile_matrix)+10)),
       axes=T, ann=T)
  apply(profile_matrix, 2, function(x) {
    lines(x, ...)
  })
}

plotProfiles <- function(profile=NULL, ci_1=NULL, ci_2=NULL, cols=NULL, xlabel=NULL,
                         ylabel=TRUE, yaxs_val=NULL, new=TRUE, smooth=FALSE, span=0.1, write=FALSE,
                         fname=NULL, main=NULL, ...) {

  if (new) {
    if (write) {
      png(file=paste(fname, ".png", sep=""))
    } else {
      x11("", 7, 7)
    }
  }
  if (is.null(yaxs_val)) {yaxs_val1 <- c((min(profile)-10), (max(profile)+10))}
  else {yaxs_val1=yaxs_val}
  plot(1, 1, type="n",
       xlim=c(1,nrow(profile)),
       ylim=yaxs_val1,
       ann=FALSE, axes=FALSE, ...)
    
  col.rgb <- col2rgb(cols, alpha=TRUE)
  col.rgb[4,] <- 50
  
  if (smooth) {
    profile <- smoothProfile(profile, span)
    if (length(ci_1) > 0) {
      ci_1 <- smoothCI(ci_1, span)
      ci_2 <- smoothCI(ci_2, span)
    }
  }
  
  for(i in 1:ncol(profile)) {
    if (length(ci_1) > 0) {
      polygon(c(c(1:nrow(profile)), c(nrow(profile):1)),
              c(ci_1[,i],
                rev(ci_2[,i])),
              col=rgb(col.rgb[1,i], col.rgb[2,i], col.rgb[3,i], alpha=col.rgb[4,i], maxColorValue=255),
              border=F)
    }
    lines(profile[,i], col=cols[i], ...)
  }
  mtext(main, main=T)
  #axis(1, at=c(0, nrow(profile)/2, nrow(profile)), labels=xlabel)
  #axis(2)
  box()
  if (ylabel) {mtext("RMS", 2, line=3, cex=1.6)}
  if (write) {dev.off()}
}

# 'sample'

plotProfiles.All <- function(AnnoSet=NULL, sample=NULL, cols=NULL, write=FALSE, ...) {
  data_path <- paste("~/analysis/profiles/", AnnoSet , sep="")
  dir_files <- list.files(data_path)
  select_files <- dir_files[grep(sample, dir_files)]
  val_files <- sort(select_files[grep("mean$", select_files)])
  ci_files_1 <- sort(select_files[grep("bootCI_1", select_files)])
  ci_files_2 <- sort(select_files[grep("bootCI_2", select_files)])
  
  val_data <- lapply(val_files, function(x) {
    d <- read.delim(paste(data_path, x, sep="/"), header=F)
    return(d)
  })
  ci_data_1 <- lapply(ci_files_1, function(x) {
    d <- read.delim(paste(data_path, x, sep="/"), header=F, colClasses="numeric")
    return(d)
  })
  ci_data_2 <- lapply(ci_files_2, function(x) {
    d <- read.delim(paste(data_path, x, sep="/"), header=F, colClasses="numeric")
    return(d)
  })
  
  for(i in 1:length(val_data)) {
    profile_name <- val_files[i]
    print(profile_name)
    fname = paste(data_path, "plots", val_files[i], sep="/")
    plotProfiles(profile=val_data[[i]], ci_1=ci_data_1[[i]],
                      ci_2=ci_data_2[[i]], cols=cols, write=write, fname=fname,
                 main=profile_name, ...)
  }
}

plotProfilesGroup <- function(profile=NULL, ci_1=NULL, ci_2=NULL, cols=NULL, xlabel=NULL, ylabel=TRUE,
                              new=TRUE, smooth=FALSE, span=0.1, write=FALSE, fname=NULL, ...) {
  if (new) {
    if (write) {
      png(file=paste(fname, ".png", sep=""))
    } else {
      x11("", 7, 7)
    }
  }  
  plot(1, 1, type="n",
       xlim=c(1,nrow(profile)),
       ylim=c((min(profile)-10), (max(profile)+10)),
       ann=FALSE, axes=FALSE, ...)
    
  col.rgb <- col2rgb(cols, alpha=TRUE)
  col.rgb[4,] <- 50
  
  if (smooth) {
    profile <- smoothProfile(profile, span)
    if (length(ci_1) > 0) {
      ci_1 <- smoothCI(ci_1, span)
      ci_2 <- smoothCI(ci_2, span)
    }
  }
  
  for(i in 1:ncol(profile)) {
    if (length(ci_1) > 0) {
      polygon(c(c(1:nrow(profile)), c(nrow(profile):1)),
              c(ci_1[,i],
                rev(ci_2[,i])),
              col=rgb(col.rgb[1,i], col.rgb[2,i], col.rgb[3,i], alpha=col.rgb[4,i], maxColorValue=255),
              border=F)
    }
    lines(profile[,i], col=cols[i], ...)
  }

  box()
  if (ylabel) {mtext("RMS", 2, line=2, cex=1.2)}
  if (write) {dev.off()}
}

plotProfilesGroup.All <- function(AnnoSet=NULL, sample=NULL, cols=NULL, write=FALSE, ...) {
  data_path <- paste("~/analysis/profiles/", AnnoSet , sep="")
  dir_files <- list.files(data_path)
  select_files <- dir_files[grep(sample, dir_files)]
  val_files <- sort(select_files[grep("mean", select_files)])
  ci_files_1 <- sort(select_files[grep("bootCI_1", select_files)])
  ci_files_2 <- sort(select_files[grep("bootCI_2", select_files)])
  
  val_data <- lapply(val_files, function(x) {
    d <- read.delim(paste(data_path, x, sep="/"), header=F)
    return(d)
  })
  ci_data_1 <- lapply(ci_files_1, function(x) {
    d <- read.delim(paste(data_path, x, sep="/"), header=F, colClasses="numeric")
    return(d)
  })
  ci_data_2 <- lapply(ci_files_2, function(x) {
    d <- read.delim(paste(data_path, x, sep="/"), header=F, colClasses="numeric")
    return(d)
  })
  
  for(i in 1:length(val_data)) {
    fname = paste(data_path, "plots", val_files[i], sep="/")
    plotProfilesGroup(profile=val_data[[i]], ci_1=ci_data_1[[i]],
                      ci_2=ci_data_2[[i]], cols=cols, write=write, fname=fname, ...)
  }
}

smoothProfile <- function(mat, span = 0.1) {
  out_mat <- apply(mat, 2, function(x) {
    out <- predict(loess(x~c(1:length(x)), span = span))
    return(out)
 })}

smoothCI <- function(mat, span = 0.1) {
  out_mat <- apply(mat, 2, function(x) {
    outs <- foreach(y=c(1,2),.combine=".wrap") %do% {
      out_vals <- .unwrap(x, y)
      out_smooth <- predict(loess(out_vals~c(1:length(out_vals)), span=span))
      return(out_smooth)
    }
    return(outs)
  })
}

.unwrap <- function(list_vec, ind) {
  vals <- unlist(lapply(list_vec, function(x) x[ind]))
  return(vals)
}

.wrap <- function(vec1, vec2) {
  out <- list()
  for(i in 1:length(vec1)) {
    out <- c(out, list(vec1[i], vec2[i]))
  }
  return(out)
}

xtss <- c("-20 kb", "TSS", "+20 kb")
xtes <- c("-20 kb", "TES", "+20 kb")
xgene <- c("-10 kb", "TSS", "TES", "+10 kb")
xmk4 <- c("-2 kb", "", "", "+2 kb")
xmk4.2 <- c("-10 kb", "", "", "+2kb")

plotTSSandTES.both <- function(annoset1, annoset2, sample, group2, fname=NULL, cols=NULL, main_lab=NULL, ...) {
  data_path <- "~/analysis/profiles/plots/"
  if (is.null(fname)) {
    x11("",10,10)
  } else {
    png(paste(data_path, fname, ".png", sep=""), width=720, height=720)
  }
  par(mfcol=c(2,2), mar=c(2, 4, 4, 0)+.1, oma=c(1,7,1,1))
  plotHMCandMC(annoset1, sample, group2, cols=cols, lab="TSS", combine=TRUE)
  plotHMCandMC(annoset2, sample, group2, cols=cols, lab="TES", ylabel=FALSE, combine=TRUE)
  par(las=1)
  mtext("5hmC", at=.73, side=2, line=0, outer=TRUE, cex=1.6, font=2)
  mtext("5mC", at=.23, side=2, line=0, outer=TRUE, cex=1.6, font=2)
  par(las=0)
  mtext(main_lab, main=TRUE, outer=TRUE, line=-2, cex=2, font=2)
  if (!is.null(fname)) {
    dev.off()
  }
}

plotHMCandMC <- function(annoset, sample, group2, fname=NULL, cols=NULL, lab=NULL,
                         ylabel=TRUE, combine=FALSE) {
  data_path <- paste("~/analysis/profiles/", annoset, sep="")
  dir_files <- list.files(data_path)
  sample_files <- dir_files[grep(group2, dir_files)]
  sample_files <- sample_files[grep(paste("^", sample, sep=""), sample_files)]
  hmc_files <- sort(sample_files[grep("_hmedip", sample_files)])
  mc_files <- sort(sample_files[grep("_medip", sample_files)])
  hmc_data <- lapply(hmc_files, function(x) {
    d <- read.delim(paste(data_path, x, sep="/"), header=F)
    return(d)
  })
  mc_data <- lapply(mc_files, function(x) {
    d <- read.delim(paste(data_path, x, sep="/"), header=F, colClasses="numeric")
    return(d)
  })
  if (!combine) {
    if (is.null(fname)) {
      x11("",10,5)
    } else {
      png(paste(fname, ".png", sep=""))
    }
    par(mfrow=c(1,2))
  }  
  .profile <- function(data) {
    plotProfilesGroup(data[[1]], data[[2]], data[[3]], cols=cols,
                    ylabel=ylabel, new=FALSE)
    l <- unlist(nrow(data[[1]]))
    a <- 1
    b <- l / 2 + 1
    c <- l
    lab_dist <- l / 2 * 200 / 1000 
    axis(1, at=c(a, b, c), labels=c(paste("+", lab_dist, " kb", sep=""),
                             lab,
                             paste("-", lab_dist, " kb", sep="")),
         cex.axis=1.2)
    axis(2, cex.axis=1.2)
  }

  .profile(hmc_data)
  .profile(mc_data)
  
  if (!is.null(fname)) {
    dev.off()
  }
}

plotTSSandTES <- function(annoset, sample, group2, fname=NULL, cols=NULL, combine=FALSE, ...) {
  data_path <- paste("~/analysis/profiles/", annoset, sep="")
  dir_files <- list.files(data_path1)
  sample_files <- dir_files[grep(group2, dir_files)]
  sample_files <- sample_files[grep(paste("^", sample, sep=""), sample_files)]
  hmc_files <- sort(sample_files[grep("_hmedip", sample_files)])
  mc_files <- sort(sample_files[grep("_medip", sample_files)])
  hmc_data <- lapply(hmc_files, function(x) {
    d <- read.delim(paste(data_path, x, sep="/"), header=F)
    return(d)
  })
  mc_data <- lapply(mc_files, function(x) {
    d <- read.delim(paste(data_path, x, sep="/"), header=F, colClasses="numeric")
    return(d)
  })
  if (is.null(fname)) {  
    if (!combine) {x11("",10,5)}
  } else {
    png(paste(fname, ".png", sep=""))
  }
  if (!combine) {par(mfrow=c(1,2))}
  
  plotProfilesGroup(tss[[1]], tss[[2]], tss[[3]], cols=cols, new=FALSE, ...)
  plotProfilesGroup(tes[[1]], tes[[2]], tes[[3]], cols=cols, new=FALSE, ...)
  if (!is.null(fname)) {
    dev.off()
  }
}

plotGenes <- function(annoset, sample, group2, fname=NULL, cols=NULL, ...) {
  data_path <- paste("~/analysis/profiles/", annoset, sep="")
  dir_files <- list.files(data_path)
  sample_files <- dir_files[grep(group2, dir_files)]
  sample_files <- sample_files[grep(paste("^", sample, sep=""), sample_files)]
  hmc_files <- sort(sample_files[grep("_hmedip", sample_files)])
  mc_files <- sort(sample_files[grep("_medip", sample_files)])
  hmc_data <- lapply(hmc_files, function(x) {
    d <- read.delim(paste(data_path, x, sep="/"), header=F)
    return(d)
  })
  mc_data <- lapply(mc_files, function(x) {
    d <- read.delim(paste(data_path, x, sep="/"), header=F, colClasses="numeric")
    return(d)
  })
  
  out_path <- paste("~/analysis/profiles/plots/", sep="")
  if (is.null(fname)) {
    x11("",9,9)
  } else {
    png(paste(out_path, fname, ".png", sep=""), width=720, height=720)
  }
  par(mfcol=c(2,1), mar=c(2, 4, 4, 0)+.1, oma=c(1,5,1,1))
  plotProfilesGroup(hmc_data[[1]], hmc_data[[2]], hmc_data[[3]], cols=cols, 
                    new=FALSE, ...)
  axis(1, at=c(0,51,100,150), labels=xgene)
  axis(2)
  par(las=1)
  mtext("5hmC", at=.73, side=2, line=0, outer=TRUE, cex=1.2, font=2)
  par(las=0)
  plotProfilesGroup(mc_data[[1]], mc_data[[2]], mc_data[[3]], cols=cols, new=FALSE, ...)
  axis(1, at=c(0,51,100,150), labels=xgene)
  axis(2)
  par(las=1)
  mtext("5mC", at=.23, side=2, line=0, outer=TRUE, cex=1.2, font=2)
  if (!is.null(fname)) {
    dev.off()
  }  
}

plotAll <- function(annoset, group2, cols) {
  cells <- c("omp","ngn","icam")
  groups <- paste(c("omp", "ngn"), group2, sep="_")
  if(annoset=="ends") {
    for(cell in cells) {
      for(group in groups) {
        fname <- paste("tss_tes", cell, "hmedip_medip", group, "group2", sep="_")
        plotTSSandTES.both("transcSS_W200N100", "transcES_W200N100", cell, group, fname, cols)
      }
    }    
  }
}
  
plotPeaks <- function(annoset, sample, fname=NULL, cols=NULL, xlabel=NULL,
                      ylabel=TRUE, ylim=NULL, combine=FALSE, ...) {
  data_path <- paste("~/analysis/profiles/", annoset, sep="")
  dir_files <- list.files(data_path)
  #sample_files <- dir_files[grep(group2, dir_files)]
  sample_files <- dir_files[grep(paste("^", sample, sep=""), dir_files)]
  hmc_files <- sort(sample_files[grep("_hmedip", sample_files)])
  mc_files <- sort(sample_files[grep("_medip", sample_files)])
  hmc_data <- lapply(hmc_files, function(x) {
    d <- read.delim(paste(data_path, x, sep="/"), header=F)
    return(d)
  })
  mc_data <- lapply(mc_files, function(x) {
    d <- read.delim(paste(data_path, x, sep="/"), header=F, colClasses="numeric")
    return(d)
  })
  #print(dir_files)
  if (!combine) {
    if (is.null(fname)) {
      x11("",10,5)
    } else {
      png(paste(fname, ".png", sep=""))
    }
    par(mfrow=c(1,2))
  }  

  plotProfiles(hmc_data[[1]], hmc_data[[2]], hmc_data[[3]], cols=cols,
                    ylabel=ylabel, new=FALSE, ylim=ylim, ...)
  axis(1, at=c(1, 11, 16, 25), labels=xmk4, cex.axis=1.6)
  axis(2, cex.axis=1.6)
  abline(v=11, lty=2, col="grey")
  abline(v=16, lty=1.8, col="grey")
  plotProfiles(mc_data[[1]], mc_data[[2]], mc_data[[3]], cols=cols,
                    ylabel=ylabel, ylim=ylim, new=FALSE, ...)
  axis(1, at=c(1, 11, 16, 25), labels=xmk4, cex.axis=1.6)
  axis(2, cex.axis=1.6)
  abline(v=11, lty=2, col="grey")
  abline(v=16, lty=2, col="grey")
  if (!is.null(fname)) {
    dev.off()
  }
}

  
plotPeaks.cells <- function(annoset, sample, fname=NULL, cols=NULL, xlabel=NULL,
                      ylabel=TRUE, yaxs_type="free", combine=FALSE, ...) {
  data_path <- paste("~/analysis/profiles/", annoset, sep="")
  dir_files <- list.files(data_path)
  sample_files <- dir_files[grep(sample, dir_files)]
  
  omp_files <- sort(sample_files[grep("omp", sample_files)])
  ngn_files <- sort(sample_files[grep("ngn", sample_files)])
  icam_files <- sort(sample_files[grep("icam", sample_files)])
  omp_data <- lapply(omp_files, function(x) {
    d <- read.delim(paste(data_path, x, sep="/"), header=F)
    return(d)
  })
  ngn_data <- lapply(ngn_files, function(x) {
    d <- read.delim(paste(data_path, x, sep="/"), header=F, colClasses="numeric")
    return(d)
  })
  icam_data <- lapply(icam_files, function(x) {
    d <- read.delim(paste(data_path, x, sep="/"), header=F, colClasses="numeric")
    return(d)
  })
  #print(dir_files)
  if (!combine) {
    if (is.null(fname)) {
      x11("",10,5)
    } else {
      png(paste(fname, ".png", sep=""))
    }
    par(mfrow=c(1,3))
  }
  yaxs_val <- NULL
  if (yaxs_type == "free") {yaxs_val <- NULL}
  if (yaxs_type == "standard") {  
    max_val <- max(c(unlist(omp_data[[1]]),
                     unlist(ngn_data[[1]]),
                     unlist(icam_data[[1]]))) + 10
    min_val <- min(c(unlist(omp_data[[1]]),
                     unlist(ngn_data[[1]]),
                     unlist(icam_data[[1]]))) - 10
    yaxs_val <- c(min_val, max_val)
  }
  .profile <- function(data, col, ylabel) {
     plotProfiles(data[[1]], data[[2]], data[[3]], cols=col,
                    ylabel=ylabel, new=FALSE, yaxs_val=yaxs_val, ...)
     l <- length(unlist(data[[1]]))
     a <- 1
     b <- (l - 5) / 2 + 1
     c <- b + 5
     d <- l
     
     axis(1, at=c(a,b,c,d), labels=xmk4.2, cex.axis=1.6)
     axis(2, cex.axis=1.6)
     abline(v=b, lty=2, col="grey")
     abline(v=c, lty=2, col="grey")
  }
  .profile(icam_data, cols[1], ylabel=TRUE)
  .profile(ngn_data, cols[2], ylabel=FALSE)
  .profile(omp_data, cols[3],ylabel=FALSE)
  if (!is.null(fname)) {
    dev.off()
  }
}

plotPeaks.all <- function(annoset, fname=NULL, cols=NULL, yaxs_type="free", ...) {
  data_path <- "~/analysis/profiles/plots/"
  if (is.null(fname)) {
    x11("",12,9)
  } else {
    png(paste(data_path, fname, ".png", sep=""), width=1080, height=720)
  }
  par(mfrow=c(2,3), mar=c(2, 4, 4, 0)+.1, oma=c(1,8,1,1))
  plotPeaks.cells(annoset, "_hmedip", cols=cols, xlabel=xmk4,
                  yaxs_type="free", combine=TRUE, ...)
  plotPeaks.cells(annoset, "_medip", cols=cols, xlabel=xmk4,
                  yaxs_type=yaxs_type, combine=TRUE, ...)
  #plotPeaks(annoset, "omp", cols=cols[3], xlabel=xmk4, ylabel=FALSE, combine=TRUE, ...)
  par(las=1)
  mtext("5hmC", at=.74, side=2, line=2, outer=TRUE, cex=1.6, font=2)
  mtext("5mC", at=.24, side=2, line=2, outer=TRUE, cex=1.6, font=2)
  par(las=0)
  mtext("ICAM", at=.19, side=3, line=-2, outer=TRUE, cex=1.6, font=2)
  mtext("Neurog1", at=.53, side=3, line=-2, outer=TRUE, cex=1.6, font=2)
  mtext("OMP", at=.86, side=3, line=-2, outer=TRUE, cex=1.6, font=2)
  if (!is.null(fname)) {
    dev.off()
  }
}

boot.samplemean <- function(x, d) {
  return(mean(x[d]))
}

AnnoSetSimp.bootCI <- function(profile_matrix, FUN=boot.samplemean) {
  require(boot)  
  trimmean <- function(x, d, trim=0) {
    return(mean(x[d], trim))
  }
  bootObj <- list(boot(data=profile_matrix[,1], statistic=FUN, R=1000))
  for (i in 2:ncol(profile_matrix)) {
    bootObj <- c(bootObj, list(boot(data=profile_matrix[,i], statistic=FUN, R=1000)))
  }
  bootCI <- lapply(bootObj, function(x) boot.ci(x, type="perc"))
  percentCI <- lapply(bootCI, function(x) x$percent[1,4:5])
  return(percentCI)
}

AnnoSetSimp.bootCImean <- function(vals) { 
 bootObj <- boot(data=vals, statistic=boot.samplemean, R=100)
  bootCI <- boot.ci(bootObj, type="perc")
  percentCI <- as.array(bootCI$percent[1,4:5])
  return(percentCI)
}

generateQuantilesAll <- function(data=NULL, probs=c(0,.25,.5,.75,1), fname=NULL) {
  it <- iter(data, by="column") 
  quantiles <- foreach(curr_it = it, .combine="cbind") %do% {
    q <- quantile(curr_it, probs)
    #return(q)
    out_mat <- rep(0, times=nrow(curr_it))
    for(i in 1:length(q-1)) {
      out_mat[curr_it > q[i] & curr_it <= q[i+1]] <- i
    }
    row.names(out_mat) <- rownames(data)
    write.table(out_mat, file=paste("~/analysis/rna/quantiles/", fname, sep=""), quote=FALSE,
                col.names=FALSE)
    return(out_mat)
  }
 
  row.names(quantiles) <- rownames(data)
  return(quantiles)  
}


generateQuantiles <- function(data=NULL, probs=c(0,.25,.5,.75,1), fname=NULL) {
    vals <- data[,2]
    q <- quantile(vals, probs)
    #return(q)
    out_mat <- rep(NA, times=nrow(data))
    out_mat[vals >= q[1] & vals <= q[2]] <- 1
    for(i in 2:length(q-1)) {
      out_mat[vals > q[i] & vals <= q[i+1]] <- i
    }
    out_mat <- data.frame(data[,1], out_mat)
    write.table(out_mat, file=paste("~/analysis/rna/quantiles/", fname, sep=""), quote=FALSE,
                sep="\t",row.names=FALSE, col.names=FALSE)
    return(out_mat)
}

collectValuesOverRegion <- function(data=NULL, reads=NULL, region=NULL, FUN=mean) {
  group <- data@set_group
  names <- data@set_names
  data_name <- unlist(strsplit(data@name, "/"))
  data_name <- data_name[length(data_name)]
  #out_path <- paste("~/analysis/profiles", data_name, sep="/")
  #if (!file.exists(out_path)) {dir.create(out_path)}
  reads <- .extract(data@read_data, reads)
  reads_vals <- getDataByType(reads, "norm")
  
  r_split_by_name <- split(reads_vals, names)

  vals <- foreach(r=r_split_by_name) %do% {
     
     r_vals <- r[region[1]:region[2]]

     vals <- do.call(FUN, list(r_vals))
     return(vals)
  }
  names(vals) <- unique(names)
  return(unlist(vals))
}

threshByCollectedVals <- function(data=NULL, reads=NULL, region=NULL, FUN=mean,
                                  thresh=0, greater=TRUE, write=FALSE, fname=NULL) {

  vals <- collectValuesOverRegion(data=data, reads=reads, region=region, FUN=mean)
  thresh_vals=NULL
  if (greater) {thresh_vals <- vals[vals >= thresh]}
  else {thresh_vals <- vals[vals <= thresh]}
  if (write) {
    anno_name <- .nameSplit(data)
    out_path <- paste("~/analysis/profiles", anno_name, "thresh", sep="/")
    if(!file.exists(out_path)) {dir.create(out_path)}
    write.table(thresh_vals, file=paste(out_path, fname, sep="/"), quote=FALSE, sep="\t",
                col.names=FALSE)
  }
  return(thresh_vals)
  
}

threshByCollectedVals.All <- function(data=NULL, region=NULL, FUN=mean,
                                      thresh_vals=c(400, 600), write=TRUE) {
  read_names <- .getNames(data@read_data)
  x <- foreach(curr_read_data=read_names) %do% {
    print(curr_read_data)
    data_name <- unlist(strsplit(curr_read_data, ".bed"))[1]
    fname <- paste(data_name, "region", region[1], region[2],
                   "below", thresh_vals[1], sep="_")
    threshByCollectedVals(data=data, reads=curr_read_data, region=region, FUN=FUN,
                                      thresh=thresh_vals[1], greater=FALSE, write=write,
                                      fname=fname)
    fname <- paste(data_name, "region", region[1], region[2],
                   "above", thresh_vals[2], sep="_")
    threshByCollectedVals(data=data, reads=curr_read_data, region=region, FUN=FUN,
                                      thresh=thresh_vals[2], greater=TRUE, write=write,
                                      fname=fname)
  } 
}

importThreshSetExportBED <- function(file, ref_annoset=NULL) {
  file_path = paste("~/analysis/profiles", .nameSplit(ref_annoset), "thresh", sep="/")
  subset <- read.delim(paste(file_path, file, sep="/"), header=FALSE)
  anno_name <- ref_annoset@set_names
  index <- anno_name%in%as.character(subset[,1])
  #return(index)
  subset_chr <- ref_annoset@set_chr[index]
  subset_start <- ref_annoset@set_pos[index]
  subset_end <- subset_start + ref_annoset@bin_size + - 1
  subset_names <- ref_annoset@set_names[index]
  subset_group <- ref_annoset@set_group[index]
  subset_strand <- ref_annoset@set_strand[index]
  
  out <- cbind(subset_chr, subset_start, subset_end,
               subset_names, subset_group, subset_strand)
  out_path <- paste(file_path, "subset", sep="")
  if(!file.exists(out_path)) {dir.create(out_path)}
  write.table(out, file=paste(out_path, paste(file, "subset.bed", sep="_"), sep="/"),
              quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
  #return(out)
}

importThreshSetExportBED.All <- function(select_string=NULL, ref_annoset=NULL) {
  file_path <- paste("~/analysis/profiles", .nameSplit(ref_annoset), "thresh", sep="/")
  dir_files <- list.files(file_path)
  
  select_files <- .selectFiles(dir_files, select_string)

  for (file in select_files) {
    importThreshSetExportBED(file, ref_annoset)
  }
  
}
