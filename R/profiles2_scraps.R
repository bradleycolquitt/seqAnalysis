### MP all MEDIPS for all annotations
### NOT FUNC
MP.all <- function(bin.size=200) {
  medips.path <- "~/data/medips"
  dir_files <- list.files(medips.path)
  select_files <- dir_files[grep(as.character(bin.size), dir_files)]
  for(file in select_files) {
    load(paste(medips.path, file, sep="/"))
    
  }
}

MP.plotHMCandMC <- function(annotation, sample, group2=NULL, fname=NULL,
                          ylabel=TRUE, stack=FALSE,
                         combine=FALSE, standard=TRUE, ...) {

  ## Define files to be plotted
  data.path <- paste("~/storage/analysis/mprofiles/", sep="")
  dir_files <- list.files(data.path)
  dir_files <- dir_files[-grep("moe", dir_files)]
  files.omp <- sort(dir_files[grep("omp", dir_files)])
  files.ngn <- sort(dir_files[grep("ngn", dir_files)])
  files.icam <- sort(dir_files[grep("icam", dir_files)])
   
  ## Get files
  data.omp <- lapply(files.omp, function(x) {
    file.path <- paste(data.path, x, "profiles", sep="/")
    profileRead(file.path, annotation, group2)
  })  
  data.ngn <- lapply(files.ngn, function(x) {
    file.path <- paste(data.path, x, "profiles", sep="/")
    profileRead(file.path, annotation, group2)
  })
  data.icam <- lapply(files.icam, function(x) {
    file.path <- paste(data.path, x, "profiles", sep="/")
    profileRead(file.path, annotation, group2)
  })
  
  ## Set up graphics device
  if (!combine) {
    if (!stack) {
      if (is.null(fname)) {
        x11("",20,10)
      } else {
        png(paste(plot.path, fname, sep="/"))
      }
      par(mfrow=c(3,2))
    } else {
      if (is.null(fname)) {
        #x11("", 12, 7)
        x11("", 12, 6)
      } else {
        #pdf(paste(plot.path, paste(fname, ".pdf", sep=""), sep="/"), width=9, height=6)
        pdf(paste(plot.path, paste(fname, ".pdf", sep=""), sep="/"), width=10, height=5)
      }
      par(mfrow=c(1,3), mar=c(3, 4, 2, 1)+.1, oma=c(1,4,2,1))
    }
  }
  
  ## Compute y axis scale, and send data to plotter
  .profile <- function(data, combine, stack, y.val, cols) {
    MP.plotAnno(data, annotation, group2, cols=col2, lab=c("",""),
                    ylabel=ylabel, combine=TRUE, y.val=y.val, stack=stack,
                standard=standard) 
  }
  
  y.val <- NULL
  if (!stack) {
    if(standard) y.val <- getRange(data.icam, buffer=50)
    sapply(c(1:2), function(x) .profile(data.icam[[x]], cols=cols[1], y.val=y.val,
                                      combine=combine, stack=stack))
    if(standard) y.val <- getRange(data.ngn, buffer=50)
    sapply(c(1:2), function(x) .profile(data.ngn[[x]], cols=cols[2], y.val=y.val,
                                      combine=combine, stack=stack))
    if(standard) y.val <- getRange(data.omp, buffer=50)
    sapply(c(1:2), function(x) .profile(data.omp[[x]], cols=cols[3], y.val=y.val,
                                      combine=combine, stack=stack))
  } else { 
    if(standard) y.val <- getRange(c(data.icam, data.ngn, data.omp), buffer=40)
    .profile(data.icam, y.val=y.val, combine=combine, stack=stack)
    #if(standard) y.val <- getRange(data.ngn)
    .profile(data.ngn, y.val=y.val, combine=combine, stack=stack)
    #if(standard) y.val <- getRange(data.omp, buffer=50)
    .profile(data.omp, y.val=y.val, combine=combine, stack=stack)
  }
  par(las=1)
  mtext("HBC", at=.19, side=3, outer=T, cex=1.6)
  mtext("GBC", at=.525, side=3, outer=T, cex=1.6)
  mtext("OSN", at=.85, side=3, outer=T, cex=1.6)
  par(las=0)
  mtext("AMS", side=2, cex=1.6, line=1, outer=T)
  if (!is.null(fname)) {
    dev.off()
  }

}

MP.plotGenes <-  function(annotation, sample, group2=NULL, fname=NULL, cols=NULL,
                          lab=NULL, standard=TRUE, combine=FALSE) {

  data.hmc <- MP.getData(mp.path, c(paste("^", sample, sep=""), "_hmedip"), "moe", annotation, group2)
  data.mc <- MP.getData(mp.path, c(paste("^", sample, sep=""), "_medip"), "moe", annotation, group2)
  
  
  data.hmc <- data.hmc[[1]]
  data.mc <- data.mc[[1]]
  data.hmc <- splitReform(data.hmc)
  data.mc <- splitReform(data.mc)
  
  if (!combine) {
    if (is.null(fname)) {
      x11("",12,6)
    } else {
      fpath <- "~/storage/analysis/mprofiles/plots/"
      pdf(paste(fpath, fname, ".pdf", sep=""), width=10, height=5)
    }
    #par(mfrow=c(2,1), mar=c(2,4,1,1)+.1, oma=c(1, 5, 1, 1))
    par(mfrow=c(1,2), mar=c(2,4,1,1)+.1, oma=c(1, 2, 4, 1))
  }  
  .profile <- function(data, y.val) {
    MP.plotAnno(data, annotation, cols=cols, lab=c("TSS", "TES"),
                    ylabel=ylabel, y.val=y.val, standard=TRUE, stack=TRUE, combine=TRUE) 
  }
  y.val <- NULL
  if(standard) y.val <- getRange(data.hmc, buffer=100)
  .profile(data.hmc, y.val)
  if(standard) y.val <- getRange(data.mc, buffer=100)
  .profile(data.mc, y.val)
  #mtext("AMS", at=.775, side=2, outer=T, line=-1, cex=1.6)
  #mtext("AMS", at=.275, side=2, outer=T, line=-1, cex=1.6)
  mtext("AMS", at=.5, side=2, outer=T, line=-3, cex=1.6)
  par(las=1)
  #mtext("5hmC", at=.775, side=2, outer=T, line=1, cex=1.6)
  #mtext("5mC", at=.275, side=2, outer=T, line=1, cex=1.6)
  mtext("5hmC", at=.275, side=3, outer=T, line=0, cex=1.6)
  mtext("5mC", at=.775, side=3, outer=T, line=0, cex=1.6)
  #mtext(paste(sample, " -> ", group2, sep=""), side=3, outer=T, cex=1.6)
  if (!is.null(fname)) {
    dev.off()
  }
}

MP.plotCells <- function(annotation, sample="", group2=NULL, fname=NULL, cols=col3,
                         lab=NULL, ylabel=TRUE, xrange=NULL, stack=TRUE,
                         standard=TRUE) {
  ## Define files to be plotted
  data.hmc <- MP.getData(mp.path, c(paste("^", sample, sep=""), "_hmc"), annotation=annotation, group2=group2)
  data.mc <- MP.getData(mp.path, c(paste("^", sample, sep=""), "_mc"), "", annotation=annotation, group2=group2)
  
  ## Set up graphics device
  
  if (!stack) {
    if (is.null(fname)) {
      x11("",20,10)
    } else {
      pdf(paste(plot.path, paste(fname, ".pdf", sep=""), sep="/"), width=8, height=6)
    }
    par(mfrow=c(2,3))
  } else {
    if (is.null(fname)) {
      x11("", 7, 10)
    } else {
      pdf(paste(plot.path, paste(fname, ".pdf", sep=""), sep="/"), width=6, height=8)
    }
    par(mfrow=c(2,1), mar=c(3, 4, 2, 1)+.1, oma=c(1,5,1,1))
  }
    
  ## Compute y axis scale, and send data to plotter
  .profile <- function(data, cols, combine, stack, y.val) {
    MP.plotAnno(data, annotation, cols=cols, lab=lab,
                y.val=y.val, combine=FALSE, stack=stack)
  }

  if (!is.null(xrange)) {
    data.hmc <- trimData(data.hmc, xrange)
    data.mc <- trimData(data.mc, xrange)
  }
  
  y.val <- NULL
  if (!stack) {
    if(standard) y.val <- getRange(data.hmc, buffer=25)
    if (length(data.hmc) > 1) {
      sapply(c(1:length(data.hmc)), function(x) .profile(data.hmc[[x]], cols=cols[x], y.val=y.val,
                                      combine=combine, stack=stack))
    } else {
      .profile(data.hmc, cols=cols, y.val=y.val, combine=combine, stack=stack)
    }
    if(standard) y.val <- getRange(data.mc)
    sapply(c(1:3), function(x) .profile(data.mc, cols=cols[x], y.val=y.val,
                                      combine=combine, stack=stack))
  } else {
    if(standard) y.val <- getRange(data.hmc, buffer=25)
    .profile(data.hmc, cols=cols, y.val=y.val, combine=combine, stack=stack)
    if(standard) y.val <- getRange(data.mc, buffer=25)
    .profile(data.mc, cols=cols, y.val=y.val, combine=combine, stack=stack)
  }
  mtext("AMS", at=.775, side=2, outer=T, line=-1, cex=1.6)
  mtext("AMS", at=.275, side=2, outer=T, line=-1, cex=1.6)
  par(las=1)
  mtext("5hmC", at=.775, side=2, outer=T, line=1, cex=1.6)
  mtext("5mC", at=.275, side=2, outer=T, line=1, cex=1.6)
  
  if (!is.null(fname)) {
    dev.off()
  }

}

MP.plot <- function(annotation, sample, group2=NULL, cols=NULL, lab=NULL, type="range") {
  ## Get sample names
  data <- MP.getData(mp.path, sample, annotation=annotation, group2=group2)[[1]]
  #return(data)
  #data <- splitReform(data)
#  return(data)
  .profile <- function(data, y.val) {
    MP.plotAnno(list(data), annotation, cols=cols, lab=lab,
                y.val=y.val) 
    }
  #y.val <- NULL
  y.val <- getRange(data, buffer=0)
  print(y.val)
  .profile(data, y.val)
  
  #mtext("AMS", at=.775, side=2, outer=T, line=-1, cex=1.6)
  #mtext("AMS", at=.275, side=2, outer=T, line=-1, cex=1.6)
  mtext("Normalized read count", at=.5, side=2, outer=T, line=-3, cex=1.6)
  par(las=1)
  #mtext("5hmC", at=.775, side=2, outer=T, line=1, cex=1.6)
  #mtext("5mC", at=.275, side=2, outer=T, line=1, cex=1.6)
  #mtext("5hmC", at=.275, side=3, outer=T, line=0, cex=1.6)
  #mtext("5mC", at=.775, side=3, outer=T, line=0, cex=1.6)
  #mtext(paste(sample, " -> ", group2, sep=""), side=3, outer=T, cex=1.6)
  #if (!is.null(fname)) {
  #  dev.off()
  #}
}

### Mean and bootstrapped 95% CI ams_A by group,
###   read in profiles from file
###   remove Inf values
###   tapply by ams_A and group
###      for boot, send to bootCI function
###   return profile

MP.makeProfile <- function(sample, mp, group2=NULL, select=NULL, cluster=NULL, write=TRUE, fname=NULL) {

  ## Read methylprofile
  file.path <- paste("~/storage/analysis/mprofiles", sample, mp, sep="/")
  mp.data <- read.delim(file.path, header=TRUE, colClasses=group.classes)
  colnames(mp.data)[c(11,20,21)] <- c("ams_A", "name", "group")

  ## Remove Inf
  mp.data <- prepMP(mp.data)

  ## Import group to split profiles
  group2.vals <- NULL
  if (!is.null(group2)) {
    group2.data <- read.delim(paste("~/analysis/rna/quantiles", group2, sep="/"),
                              header=FALSE)
    group2.data <- group2.data[group2.data[,2]>0,]
    group2.vals <- group2.data[match(mp.data$name, group2.data[,1]), 2]
  }
  
  ## Select observations
  if (!is.null(select)) {
    select.path <- paste(mp.path, "features", select, sep="/") 
    select.data <- read.delim(select.path, header=TRUE)
    select.data <- select.data[!duplicated(select.data[,1]),]
    rownames(select.data) <- select.data[,1]
    mp.data <- mp.data[mp.data$name%in%rownames(select.data),]
  }

  ## Select observations based on previous clustering
  if (!is.null(cluster)) {
    select.path <- paste(mp.path, "features/clustering/hclust", cluster, sep="/") 
    select.data <- scan(file=select.path, what=character(0))
    #select.data <- select.data[!duplicated(select.data[,1]),]
    #rownames(select.data) <- select.data[,1]
    mp.data <- mp.data[mp.data$name%in%select.data,]
  }
  
  ## calculate profile mean
  profile <- NULL
  if (is.null(group2)) {profile <- profileCompute(mp.data, list(ams_A, list(group), mean)) 
  } else {profile <- profileCompute(mp.data, list(ams_A, list(group, group2.vals), mean))}
  #if (is.null(group2)) {profile <- with(mp.data, tapply(ams_A, group, mean))
  #} else {profile <- with(mp.data, tapply(ams_A, list(group, group2.vals), mean))}

  ## calculate bootstrapped CI
  CI <- NULL
  if (is.null(group2)) {CI <- profileCompute(mp.data, list(ams_A, list(group), bootCI))
  } else {CI <- profileCompute(ams_A, list(group, group2.vals), bootCI)}                      
  #if (is.null(group2)) {CI <- with(mp.data, tapply(ams_A, group, bootCI))
  #} else {CI <- with(mp.data, tapply(ams_A, list(group, group2.vals), bootCI))}
  
  if (write) {
    out.path <- paste("~/storage/analysis/mprofiles", sample, "profiles", sep="/")
    writeModule(out.path, profile, CI, group2)
    #if (!file.exists(out.path)) {
    #  dir.create(out.path)
    #}
    #out.name <- NULL
    #if (!is.null(group2)) {out.name <- paste(mp, group2, sep="_")
    #} else if (!is.null(select)) {out.name <- paste(mp, fname, sep="_")
    #} else if (!is.null(cluster)) {out.name <- paste(mp, cluster, sep="_")}
    #else {out.name <- mp}                              
    #write.table(as.matrix(profile),
    #            file=paste(out.path, paste(out.name, "mean", sep="_"), sep="/"),
    #            quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)     
    
    #lapply(c(1:2), function(x) {
    #    unwrapped <- NULL
    #    if (is.null(group2)) {unwrapped <- .unwrap(CI, x)
    #    } else {unwrapped <- apply(CI, 2, function(y) .unwrap(y, x))}
    #    write.table(unwrapped,
    #                file=paste(out.path, paste(out.name, "mean_bootCI", x, sep="_"),
    #                  sep="/"),
    #                quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
    #})
  } else {
    return(profile)
  } 
}

## Methylprofile over windowed annotation file
### PARAM: data = MEDIPS object, roi = annotation file
### RETURN: methylprofile matrix with names and group number of annotation
MP.group <- function(data, input=NULL, roi) {
  roi.data <- read.delim(paste(anno.path, roi, sep=""), header=FALSE)
  roi.input <- cbind(roi.data[,1:3], c(1:nrow(roi.data)))
  if (!is.null(input)) {mp <- MEDIPS.methylProfiling(data1=data, input=input, ROI_file=roi.input, transf=FALSE)
  } else {mp <- MEDIPS_mod.methylProfiling(data1=data, ROI_file=roi.input, transf=FALSE)}
  mp <- cbind(mp, name=roi.data[,4], group=roi.data[,5])
  return(mp)
}


## MP.group for all annotations
MP.group.allAnno <- function(data=NULL, input=NULL) {
  #anno.path <- "~/lib/annotations"
  mp.path <- paste("~/storage/analysis/mprofiles", data@sample_name, sep="/")
  if(!file.exists(mp.path)) dir.create(mp.path)
  files <- list.files(anno.path)
  files.existing <- list.files(mp.path)
  files <- files[!files%in%files.existing]
  registerDoMC(cores=1)
  mps <- foreach (file=files) %dopar% {
    cat(paste("-- ", file, " --\n", sep=""))
    mp <- MP.group(data, input, file)
    write.table(mp, file=paste(mp.path, file, sep="/"), quote=FALSE, sep="\t",
                row.names=FALSE, col.names=TRUE)
    return(mp) 
  }
}

profileCompute <- function(data, param) {
  grouping <- lapply(param[[2]], function(x) if(!is.null(x)) data[,x])
  return(tapply(data[,param[[1]]], grouping, param[[3]]))
}


MP.group2.allSamples <- function(set="d3a", anno.path, roi_name, select=4) {
  if (set=="d3a") {
    samples <- samples.d3a
  } else if (set=="cells") {
    samples <- samples.cells
  }
  for (sample in samples) {
    files.existing <- list.files(paste(profile2.path, sample, sep="/"))
    if (!roi_name %in% files.existing) {
      print(sample)
      subsetByROI.par(sample, anno.path, roi_name, select)
    }
  }  
}

MP.group2.all <- function(set="d3a", anno_set="std", select=4) {
  if (anno_set == "std") {
    anno.path <- anno.path   
  } else if (anno_set == "hires") {
    anno.path <- anno_hires.path
  }
  rois <- list.files(anno.path)
  print(rois)
  for (roi in rois) {
    print(paste("--", roi, "--", sep="" ))
    MP.group2.allSamples(set=set, anno.path=anno.path, roi_name=roi, select=select)
  }
}

MP.update <- function(anno, dummy=0) {
  data <- read.delim(anno, colClasses=c("character", "numeric", "numeric", "character", "numeric", "character", "numeric"))
  out <- cbind(data, dummy)
  names(out) <- c(colnames(data)[1:6], "raw", "norm")
  write.table(out, anno, quote=FALSE, sep="\t", row.names=FALSE)
}

MP.update.anno <- function(sample, dummy=0) {
  files <- list.files(sample)
  for (file in files) {
    print(paste("--", file, "__", sep=""))
    MP.update(paste(sample, file, sep="/"), dummy)
  }
}
MP.update.all <- function(set="d3a", dummy=0) {
  if (set=="d3a") samples <- samples.d3a
  for (sample in samples) {
    print(sample)
    path <- paste(profile2.path, sample, sep="/")
    MP.update.anno(path, dummy=dummy)
  }
}

MP.nameChange <- function() {
  
  for (sample in samples.d3a) {
    print(sample)
    sample.path <- paste(profile2.path, sample, sep="/")
    files <- list.files(sample.path)
    ind <- grep("profiles", files)
    if (length(ind) > 0) files <- files[-ind]
    
    for (file in files) {
      print(file)
      file.path <- paste(sample.path, file, sep="/")
      data <- read.delim(file.path)
      colnames(data)[7] <- "raw"
      write.table(data, file=file.path, quote=F, sep="\t", row.names=F)
    }
  }
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

prepMP <- function(vals) {
  vals <- vals[is.finite(vals$ams_A),]
  vals <- vals[!is.na(vals$ams_A),]
  return(vals)
}

