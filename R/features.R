source("~/src/R/paths.R")
source("~/src/R/blocks.R")
source("~/src/R/profiles2.R")

feature.path <- "~/lib/features_general"


makeFeatureMatrix <- function(feature, set="cells", value_type = "raw", write=TRUE) {
  if (set=="cells") {
    samples <- samples.cells
  } else if (set=="cells_rlm") {
    samples <- samples.cells.rlm
  }
  if (value_type == "raw") {
    value_column <- 5
  } else {
    value_column <- 6
  }
 
 vals <- foreach (sample=samples, .combine="cbind") %dopar% {
    data <- read.delim(paste(feature2.path, sample, feature, sep="/"), header=FALSE)
    val <- data[, value_column]
    names(val) <- data[,4]
    return(val)
  }
  
 colnames(vals) <- samples
 if (write) write.table(vals, file=paste(feature2.path, "summaries", paste(set, feature, value_type, sep="_"), sep="/"),
                        quote=FALSE, sep="\t")
 return(vals)
}

makeFeatureMatrix.all <- function(set="cells", value_type="raw") {
  files <- list.files(feature.path)
  files <- files[grep("chr", files)]
  for(file in files) {
    print(file)
    if (file.exists(paste(feature2.path, "summaries", paste(set, file, value_type, sep="_"), sep="/"))) {
      print("File exists")
      return
    }  
    a <- makeFeatureMatrix(file, set=set, value_type=value_type, write=TRUE)
  }
}

makeFeatureMatrix2 <- function(feature, set = "all", write=TRUE) {
  if (set=="cells") {
    samples <- samples.cells_norm
  } else if (set=="d3a") {
    samples <- samples.d3a_norm
  } else if (set=="all") {
    samples <- list.files(paste(feature_norm_path, feature, sep="/")) 
  }
  print(samples)
  vals <- foreach (sample=samples, .combine="cbind") %dopar% {
    print(paste(feature_norm_path, feature, sample, sep="/"))
    data <- read.delim(paste(feature_norm_path, feature, sample, sep="/"), header=FALSE)
    val <- data[,5]
    names(val) <- data[,4]
    return(val)
  }
  colnames(vals) <- samples
  if (write) write.table(vals, file=paste(feature_norm_path, "summaries",
                                 paste(set, feature, sep="_"), sep="/"),
                        quote=FALSE, sep="\t")
  return(vals)
  
  
}
prepForHeatmap <- function(mat, var_ind=c(1:3), N=1000) {
  mat <- apply(mat, 2, pseudoCountNorm)
  mat <- na.omit(mat)
  #mat <- log(mat, 2)
  mat_extremes <- quantile(mat, probs=c(.1,.9))
  mat <- mat[as.logical(apply(mat >= mat_extremes[1] & mat <= mat_extremes[2], 1, prod)),]
  mat.var <- apply(mat, 1, function(x) var(x[var_ind]))
  mat.var.q <- quantile(mat.var, probs=(length(mat.var) - N)/ length(mat.var))
  mat <- mat[mat.var >= mat.var.q,]
  return(mat)
}
makeFeatureDF <- function(set="cells", value_type="raw", FUN=mean) {
  files <- list.files(paste(feature2.path, "summaries", sep="/"))
  files <- files[grep(set, files)]
  files <- files[grep(value_type, files)]
  files_names <- lapply(files, str_split, "_")
  print(files)
  #return(files_names)
  files_names <- lapply(files_names, function(x) {
    #print(x)
    sel <- x[[1]][-grep(set, x[[1]])]
    #return(sel)
    sel <- sel[-grep(value_type, sel)]
    sel <- sel[-grep("chr", sel)]
    return(sel)
  })
  #print(files_names)
  #return(files_names)
  files_names <- unlist(lapply(files_names, paste, collapse="_"))
  stat <- list()
  ci <- list()
  data <- foreach(file=files) %do% {
    tmp <- read.delim(paste(feature2.path, "summaries", file, sep="/"))
    stat <- c(stat, list(apply(tmp, 2, FUN)))
    ci_tmp <- apply(tmp, 2, bootCI)
    rownames(ci_tmp) <- c("lower", "upper")
    ci <- c(ci, list(ci_tmp))
  }
  #return(stat)
  #return(ci)
  #return(data)
  #names(stat) <- files_names
  #names(ci) <- files_names
  stat.df <- ldply(stat)
  ci.df <- ldply(ci)
  #return(ci.df)
  #return(stat.df)
  stat.df$value_type <- "stat"
  stat.df$feature <- files_names
  fl <- length(unique(stat.df$feature))
  print(fl)
  stat.melt <- melt(stat.df)
  stat.melt$celltype <- rep(c("mOSN", "GBC", "HBC"), each=fl)
  stat.melt$modification <- rep(c("5hmC", "5mC"), each=fl * 3)
  ci.df$value_type<- c("lower", "upper")
  ci.df$feature <- rep(files_names, each=2)
  ci.melt <- melt(ci.df)
  ci.melt$celltype <- rep(c("mOSN", "GBC", "HBC"), each=fl * 2)
  ci.melt$modification <- rep(c("5hmC", "5mC"), each=fl * 6)
  #return(ci.df)
  #return(data)
  #data.df <- ldply(data)
  #return(data.df)
  #stat.melt <- melt(stat.df)
  #return(stat.melt)
  #ci.melt <- melt(ci.df)
  data.melt <- rbind(stat.melt, ci.melt)
  #data.df <- melt(data.df)
  return(data.melt)
  
}

compareFeatures <- function(feature, sample_list, value_type="raw") {
  data <- lapply(sample_list, function(sample) {
    read.delim(paste(feature2.path, sample, feature, sep="/"), header=FALSE)
  })
  
}

.sem <- function(vals) {
  sem <- sd(vals)/sqrt(length(vals) - 1)
  return(sem)
}
