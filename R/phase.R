source("~/src/seqAnalysis/R/paths.R")

loadPhaseData <- function() {
  files <- list.files(phase.path)
  print(files)
  data <- lapply(files, function(x) scan(file=paste(phase.path, x, sep="/")))
  names(data) <- files
  return(data)
}

normPhaseData <- function(data, num_regions=NULL) {
  #omp_data <- data[grep("omp_nuc", names(data))]
  #icam_data <- data[grep("icam_nuc", names(data))]
  #omp_count <- 73284106
  #icam_count <- 106293504
  #omp_data <- lapply(omp_data, function(x) 1E6 * x/omp_count)
  #icam_data <- lapply(icam_data, function(x) 1E6 * x/icam_count)
 # data <- c(omp_data, icam_data)
  data <- lapply(data, function(x) x/sum(x))
  #data <- lapply(1:length(data), function(x) data[[x]]/num_regions[x]) 
  return(data)
}

plotPhase <- function(data, phase_ind, ...) {
  plot(phase_ind, data[phase_ind], type="l",...)
} 

plotPhase.pair <- function(data, phase_ind, fname=NULL, ...) {
  if (is.null(fname)) {
    x11()
  } else {
    pdf(file=paste("~/s2/analysis/nuc/plots", fname, sep="/"), 10, 4)
  }  
  par(mfrow=c(1,4))
  for(i in 1:4) {
    plot(phase_ind, data[[i]][phase_ind], type="l", col="dark blue", ...)
    lines(phase_ind, data[[i+4]][phase_ind], col="dark red")
  }

  if (!is.null(fname)) dev.off()
}
