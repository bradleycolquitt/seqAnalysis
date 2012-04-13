source("~/src/seqAnalysis/R/paths.R")
source("~/src/seqAnalysis/R/seqUtil.R")
registerDoMC(cores=6)

bugn <- c("dark blue", "dark green")
reds <- c("dark orange", "dark red")
blgr <- c("dimgrey", "dark green")
blred <- c("dimgrey", "dark red")

ind16 <- seq(1, 2000, 16)

loadPhaseData <- function(path, filter=NULL, neg_filter=NULL) {
  files <- list.files(path)
  if (!is.null(filter)) {
    ind <- grep(filter, files)
    if (length(ind) > 0 ) files <- files[ind]
  }
  if (!is.null(neg_filter)) {
    for(i in 1:length(neg_filter)) {
      ind <- grep(neg_filter[i], files)
      if (length(ind) > 0) files <- files[-ind]
    }  
  }
  
  print(files)
  data <- lapply(files, function(x) scan(file=paste(path, x, sep="/")))
  names(data) <- files
  return(data)
}

normPhaseData <- function(data, phase_ind, num_regions=NULL) {
  #omp_data <- data[grep("omp_nuc", names(data))]
  #icam_data <- data[grep("icam_nuc", names(data))]
  #omp_count <- 73284106
  #icam_count <- 106293504
  #omp_data <- lapply(omp_data, function(x) 1E6 * x/omp_count)
  #icam_data <- lapply(icam_data, function(x) 1E6 * x/icam_count)
 # data <- c(omp_data, icam_data)
  data <- lapply(data, function(x) x/sum(x))
  
  data <- lapply(data, function(x) {
    
    norm <- mean(x[phase_ind[(length(phase_ind) - 20):length(phase_ind)]])
    print(norm)
    return(x / norm)
  })  
  #data <- lapply(1:length(data), function(x) data[[x]]/num_regions[x]) 
  return(data)
}

phase.fft <- function(data, phase_ind) {
  pd_fft <- fft(data[phase_ind[2:(length(phase_ind) - 1)]])
  pd_mag <- Mod(pd_fft)
  pd_mag <- pd_mag / min(pd_mag)
  return(pd_mag)
}

makeFFTmatrix <- function(data_path, filter=NULL, neg_filter=NULL, sort_ind=5, phase_ind=ind16) {
  pd <- loadPhaseData(path=data_path, filter=filter, neg_filter=neg_filter)
  if (sort_ind > 0) pd <- sortByName(pd, ind_pos=sort_ind)
  pd <- normPhaseData(pd, phase_ind=phase_ind)
  pd.fft <- lapply(pd, function(x) phase.fft(x, phase_ind=phase_ind))
  pd.fft.mat <- do.call("rbind", pd.fft)
  return(pd.fft.mat)
}

phase.plot.fft <- function(fft_data, step, ...) {
  plot(fft_data[2:length(fft_data)/2], type="l", axes=T, ... )
  #print(length(fft_data)/2)
  #axis(1, at=seq(1, length(fft_data)/2, 10), labels=seq(0, length(fft_data)/2 * 10, 100))
  #axis(2)
}

phase.fftMax <- function(data, trim) {
  maxes <- lapply(data, function(x) max(x[trim[1]:trim[2]]))
  return(unlist(maxes))
}

phase.fftWhichMax <- function(data, trim) {
  maxes <- lapply(data, function(x) which.max(x[trim[1]:trim[2]]))
  return(unlist(maxes))
}

plotPhase <- function(data, phase_ind, ...) {
  plot(phase_ind, data[phase_ind], type="l",...)
} 

plotPhase.pair <- function(data, phase_ind, step, fname=NULL, cols=greens, ...) {
  
  #if (is.null(fname)) {
  #  x11()
  #} else {
  #  pdf(file=paste("~/s2/analysis/nuc/plots", fname, sep="/"), 12, 4)
  #}
  #par(mfrow=c(1,4))
  
  #for(i in 1:(length(data)/2)) {
  for (i in seq(1, length(data), 2)) {
    plot(phase_ind, data[[i]][phase_ind], type="l", col=cols[1],...)
    lines(phase_ind, data[[i+step]][phase_ind], col=cols[2])
    #box()
    #axis(1)
    #axis(2, at=ylim)
  }

  if (!is.null(fname)) dev.off()
}

plotPhase.triple <- function(data, phase_ind, step, ...) {
  for (i in 1:(length(data)/3)) {
    plot(phase_ind, data[[i]][phase_ind], type="l", col="dark red", ...)
    lines(phase_ind, data[[i + step]][phase_ind], col="dark green")
    lines(phase_ind, data[[i + 2 * step]][phase_ind], col="dark blue")
  }
}

plotPhase.four <- function(data, phase_ind, step, ...) {
  for(i in seq(1, length(data), 4)) {
    #print(names(data)[[i]])
    #print(names(data)[[i+1]])
    #print(names(data)[[i+2]])
    #print(names(data)[[i+3]])
    plot(phase_ind, data[[i]][phase_ind], type="l", col=col4_mod[1], ...)
    lines(phase_ind, data[[i+1]][phase_ind], col=col4_mod[2])
    lines(phase_ind, data[[i+2]][phase_ind], col=col4_mod[3])
    lines(phase_ind, data[[i+3]][phase_ind], col=col4_mod[4])
  }
}

#total_reads <- list(omp=73284106, icam=106293504)
total_reads <- list(omp=128485884, icam=142654159)
loadPileupData <- function(path) {
  files <- list.files(path)
  print(files)
  data <- lapply(files, function(x) read.delim(file=paste(path, x, sep="/"), sep=" ", header=FALSE))
  names(data) <- files
  return(data)

}

meanPileupData <- function(data_list) {
  #out <- foreach(x=1:length(data_list)) %dopar% {
  out <- lapply(1:length(data_list), function(x) {
    data <- data_list[[x]]
    nuc_name <- str_split(names(data_list)[[x]], "_")[[1]][1]
    print(names(data_list[x]))
    data <- 10^6 * data / total_reads[[nuc_name]]
    means <- apply(data, 2, mean)
    return(means)
  })
  names(out) <- names(data_list)
  return(out)
  #data <- 10^6 * data / 
  #means <- apply(data, 2, mean)
  #means <- means/sum(means)
  #return(means)
}

smoothPileupData <- function(data, span=.05) {
  #print(length(data))
  vals <- predict(loess(data ~ c(1:length(data)), span=span))
  return(vals)
}

fitSpline <- function(data, spar=0.05) {
  return(smooth.spline(1:length(data), data, spar=spar)$y)
  
}

plotPileup.four <- function(data, ...) {
  for(i in seq(1, length(data), 4)) {
    print(names(data)[[i]])
    print(names(data)[[i+1]])
    print(names(data)[[i+2]])
    print(names(data)[[i+3]])
    plot(data[[i]], type="l", col=col4[1], ...)
    lines(data[[i+1]], col=col4[2])
    lines(data[[i+2]], col=col4[3])
    lines(data[[i+3]], col=col4[4])
  }
}

plotPileup.two <- function(data, ...) {
  for(i in 1:(length(data)/2)) {
    print(names(data)[i])
    print(names(data)[i+2])
    plot(data[[i]], type="l", col=col2[1], ...)
    lines(data[[i+2]], col=col2[2])
  }
}
