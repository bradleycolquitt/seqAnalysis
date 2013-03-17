library(RColorBrewer)
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
normPhaseData <- function(data, phase_ind=ind16, num_regions=NULL) {
  data <- lapply(data, function(x) x/sum(x))
  data <- lapply(data, function(x) {
    norm <- mean(x[phase_ind[(length(phase_ind) - 20):length(phase_ind)]])
    print(norm)
    return(x / norm)
  })  
  return(data)
}

# Fourier transform nucleosome phasing data, return magnitudes
phase.fft <- function(data, phase_ind=ind16, start=5) {
  pd_fft <- fft(data[phase_ind[start:(length(phase_ind) - 1)]])
  pd_mag <- Mod(pd_fft)
  pd_mag <- pd_mag / min(pd_mag)
  return(pd_mag)
}

# Input phasogram output directory
# Sort_ind : field number in "_" delimited name to sort by
# Useful for constructing FPKM ordered phase data
makeFFTmatrix <- function(data_path, filter=NULL, neg_filter=NULL, sort_ind=5, phase_ind=ind16) {
  pd <- loadPhaseData(path=data_path, filter=filter, neg_filter=neg_filter)
  if (sort_ind > 0) pd <- sortByName(pd, ind_pos=sort_ind)
  pd <- normPhaseData(pd, phase_ind=phase_ind)
  pd.fft <- lapply(pd, function(x) phase.fft(x, phase_ind=phase_ind))
  pd.fft.mat <- do.call("rbind", pd.fft)
  return(pd.fft.mat)
}

phase.plot.fft <- function(fft_data, start=2, ...) {
  plot(fft_data[start:length(fft_data)/2], type="l", axes=T, ... )
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
  for(i in 1:(length(data)/2)) {
  #for (i in seq(1, length(data), 2)) {
    plot(phase_ind, data[[i]][phase_ind], type="l", col=cols[1],...)
    lines(phase_ind, data[[i+step]][phase_ind], col=cols[2])
    #box()
    #axis(1)
    #axis(2, at=ylim)
  }

  if (!is.null(fname)) dev.off()
}

plotPhase.triple <- function(data, phase_ind, step, cols, ...) {
  for (i in 1:(length(data)/3)) {
    plot(phase_ind, data[[i]][phase_ind], type="l", col=cols[i], ...)
    lines(phase_ind, data[[i + step]][phase_ind], col=cols[i+1])
    lines(phase_ind, data[[i + 2 * step]][phase_ind], col=cols[i+2])
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

# Report locations of nucleosome phasing peaks
reportPhasePeaks <- function(data, phase_ind) {
  peaks <- c()
  for (ind in 3:(length(phase_ind)-2)) {
    curr_ind <- phase_ind[(ind-2):(ind+2)]
    if (data[curr_ind[1]] < data[curr_ind[2]] & data[curr_ind[2]] < data[curr_ind[3]])
      if (data[curr_ind[3]] > data[curr_ind[4]] & data[curr_ind[4]] > data[curr_ind[5]]){
        peaks <- c(peaks, curr_ind[3])}
  }
  return(peaks)
}

# Determine slope from peak data
slopeFromPeaks <- function(peaks) {
  peak_data <- data.frame(peaks=peaks, index=1:length(peaks))
  model <- lm(peaks~index, data=peak_data)
  return(coefficients(model)[2])
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
