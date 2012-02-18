source("~/src/seqAnalysis/R/paths.R")
source("~/src/seqAnalysis/R/blocks.R")
source("~/src/seqAnalysis/R/profiles2.R")
source("~/src/seqAnalysis/R/boot.R")
source("~/src/seqAnalysis/R/modeling.R")

registerDoMC(cores=10)

feature.path <- "~/lib/features_general"
feature_norm_path <- "~/s2/analysis/features/norm"

features_toplot <- factor(1:10, labels=c("refgene_1to3kb_up_chr", "cgi_chr", "Refgene_5_UTR_chr",
                                  "Refgene_CDS_chr", "Refgene_3_UTR_chr", "Refgene_exons_chr",
                                  "Refgene_intron_chr", "omp_mk4_intergenic_inter_cons.bed_chr",
                                  "phastCons30way_intergenic_merge500_thresh500_chr",
                                  "intergenic_sub_rmsk_chr"))

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
      next
    }  
    a <- makeFeatureMatrix(file, set=set, value_type=value_type, write=TRUE)
  }
}

## Make summaries
makeFeatureMatrix2 <- function(feature, set = "all", value_type, write=TRUE) {
  if (set=="cells_norm" | set=="cells") {
    samples <- samples.cells_norm
  } else if (set=="unnorm") {
    samples <- samples.cells_raw
  } else if (set=="medips_rf_1" | set=="medips_rf_2") {
    samples <- samples.medips_rf
  } else if (set=="rpm_avg_2") {
    samples <- samples.cells_norm
  } else if (set=="tfo_unnorm") {
    #value_type <- "unnorm_2"
    samples <- samples.tfo_unnorm
  } else if (set=="tfo") {
    samples <- samples.tfo
  } else if (set=="d3a") {
    samples <- samples.d3a_norm
  } else if (set=="d3a_unnorm") {
    samples <- samples.d3a_norm
  } else if (set=="all") {
    samples <- list.files(paste(feature_norm_path, feature, sep="/")) 
  } else if (set=="d3a_mrna") {
    samples <- c("moe_d3a_wt_mrna_rpkm", "moe_d3a_ko_mrna_rpkmc")
  }
  print(samples)
  vals <- foreach (sample=samples, .combine="cbind") %dopar% {
    print(paste(feature_norm_path, value_type, feature, sample, sep="/"))
    data <- read.delim(paste(feature_norm_path, value_type, feature, sample, sep="/"), header=FALSE)
    val <- data[,5]
    names(val) <- data[,4]
    return(val)
  }
  colnames(vals) <- samples
  out_path <- paste(feature_norm_path, value_type, "summaries", sep="/")
  if (!file.exists(out_path)) dir.create(out_path)
  
  if (write) write.table(vals, file=paste(feature_norm_path, value_type, "summaries",
                                 paste(set, feature, sep="_"), sep="/"),
                        quote=FALSE, sep="\t")
  return(vals)
}

makeFeatureMatrix2.all <- function(set, value_type) {
  files <- list.files(feature.path)
  files <- files[grep("chr", files)]
  for(file in files) {
    print(file)
    if (file.exists(paste(feature_norm_path, value_type, "summaries", paste(set, file, sep="_"), sep="/"))) {
      print("File exists")
      next
    }  
    a <- makeFeatureMatrix2(file, set=set, value_type=value_type, write=TRUE)
  }
}


threshMatByQ <- function(data, q) {
  
}
statCollect <- function(set, value_type, feature, transf=NULL) {
  data <- read.delim(paste(feature_norm_path, value_type, "summaries",
                           paste(set, feature, sep="_"), sep="/"),
                     header=TRUE, row.names=NULL)
  data_melt <- melt(data)
  data_melt$feature <- feature
  labels <- sapply(colnames(data)[2:ncol(data)], str_split, "_")
  celltype <- unique(unlist(lapply(labels, function(x) x[1])))
  ip <- unique(unlist(lapply(labels, function(x) x[2])))
#  print(celltype)
#  print(ip)
  data_melt$celltype <- rep(celltype, each=nrow(data))
  data_melt$ip <- rep(ip, each=nrow(data) * 2)
  return(data_melt)
}


statSummary <- function(set, value_type, feature, transf=NULL, FUN, ...) {
  data <- read.delim(paste(feature_norm_path, value_type, "summaries",
                           paste(set, feature, sep="_"), sep="/"),
                     header=TRUE, row.names=NULL)
  if (!is.null(transf)) {
    data[,2:ncol(data)] <- apply(data[, 2:ncol(data)], 2, function(x) do.call(transf, list(x)))
  }
  return(apply(data[,2:ncol(data)], 2, FUN, ...))
}

statSummary.all <- function(set, value_type, toplot=TRUE, action="summary", transf=NULL, FUN=mean, ...) {
  if (toplot) {
    features <- as.character(features_toplot)
  } else {

  
    features <- list.files(paste(feature_norm_path, value_type, "summaries", sep="/"))
                                        #print(features)
    features <- lapply(features, str_split, paste(set, "_", sep=""))
                                        #print(features)
    features <- na.omit(unlist(lapply(features, function(x) x[[1]][2])))
                                        #print(features)
  }
  
  out <- foreach (feature=features, .combine="rbind") %dopar% {
    print(feature)
    if (action=="summary") {
      result <- statSummary(set, value_type, feature, transf=transf, FUN, ...)
    } else if (action=="collect") {
      result <- statCollect(set, value_type, feature)
    }
    return(result)
  }
  if (action=="summary") rownames(out) <- features
  return(out)
}

statSummary.allCI <- function(set, value_type, transf=NULL, boot=TRUE, ...) {
  stat <- melt(statSummary.all(set=set, value_type=value_type, transf=transf, FUN=mean, ...))
  if (boot) {
    CI <- bootCI
  } else {
    CI <- ci95
  }  
  CIup <- melt(statSummary.all(set=set, value_type=value_type, transf=transf,
                               FUN=CI, stat=boot.samplemean, bound="upper" ))
  CIlow <- melt(statSummary.all(set=set, value_type=value_type, transf=transf,
                                FUN=CI, stat=boot.samplemean, bound="lower"))
  colnames(stat) <- c("feature", "sample", "median")
  colnames(CIup) <- c("feature", "sample", "CIup")
  colnames(CIlow) <- c("feature", "sample", "CIlow")
  stat_ci <- merge(stat, CIup, all.y=TRUE)
  stat_ci <- merge(stat_ci, CIlow, all.y=TRUE)
  return(stat_ci)
}

statSummary.allSEM <- function(set, value_type, FUN=sem) {
  stat <- melt(statSummary.all(set=set, value_type=value_type, FUN=median))
  SEM <- melt(statSummary.all(set=set, value_type=value_type, FUN=FUN))
  colnames(stat) <- c("feature", "sample", "median")
  colnames(SEM) <- c("feature", "sample", "SEM")
  return(merge(stat, SEM, all.y=TRUE))
}

statSummary.allSD <- function(set, value_type, trim=.05, ...) {
  stat <- melt(statSummary.all(set=set, value_type=value_type, FUN=mean, trim=trim))
  sd <- melt(statSummary.all(set=set, value_type=value_type, FUN=sd_trim, trim=trim))
  colnames(stat) <- c("feature", "sample", "median")
  colnames(sd) <- c("feature", "sample", "sd")
  return(merge(stat, sd, all.y=TRUE))
}

statSummary.allIQR <- function(set, value_type, transf=NULL, ...) {
  stat <- melt(statSummary.all(set=set, value_type=value_type, transf=transf, FUN=median, ...))
 
  up <- melt(statSummary.all(set=set, value_type=value_type, transf=transf,
                               FUN=iqr, bound="upper" ))
  low <- melt(statSummary.all(set=set, value_type=value_type, transf=transf,
                                FUN=iqr, bound="lower"))
  colnames(stat) <- c("feature", "sample", "median")
  colnames(up) <- c("feature", "sample", "up")
  colnames(low) <- c("feature", "sample", "low")
  stat_ci <- merge(stat, up, all.y=TRUE)
  stat_ci <- merge(stat_ci, low, all.y=TRUE)
  return(stat_ci)
}

sem <- function(vals) {
  return(sd(vals)/(length(vals - 1)))
}

sd_trim <- function(vals, trim=.05) {
  q <- quantile(vals, probs=c(trim, 1-trim))
  vals <- vals[vals >= q[1] & vals <= q[2]]
  return(sd(vals))
}

ci95 <- function(vals, bound="both") {
  normfit <- fitdistr(vals, "normal")
  val_norm <- qnorm(c(.025, .975), mean=normfit$estimate[1], sd=normfit$estimate[2])
  if (bound=="both") {
    return(val_norm)
  } else if (bound=="upper") {
    return(val_norm[2])
  } else if (bound=="lower") {
    return(val_norm[1])
  }
}

iqr <- function(vals, bound="both") {
  stat <- boxplot.stats(vals)$stats[c(2,4)]
  if (bound=="both") {
    return(stat)
  } else if (bound=="upper") {
    return(stat[2])
  } else if (bound=="lower") {
    return(stat[1])
  }
}

trimOutliers <- function(mat, trim) {
  out <- apply(mat, 2, function(x) {
    q <- quantile(x, probs=c(trim, 1-trim))
    return(x[x>=q[1] & x<=q[2]])
  })  
  return(out)
}
splitSummaryByQ <- function(data, q, fname=NULL) {
  qs <- apply(data, 2, quantile, probs=q)
  
    data_down <- lapply(1:ncol(data), function(x) {
      data[data[,x]<=qs[1,x],]
    })
    data_up <- lapply(1:ncol(data), function(x) {
      data[data[,x]>=qs[2,x],]
    })
  names(data_down) <- colnames(data)
  names(data_up) <- colnames(data)
  if (!is.null(fname)) {
    for(i in 1:length(data_down)) {
      
        write.table(data_down[[i]], file=paste(fname, names(data_down)[i], q[1], sep="_"),
                    quote=FALSE, sep="\t")
        write.table(data_up[[i]], file=paste(fname, names(data_up)[i], q[2], sep="_"),
                    quote=FALSE, sep="\t")
      
    }
  }
  out <- list(data_down, data_up)
  names(out) <- q
  return(out)
}

saveBedBySubset <- function(data, bed, fname=NULL) {
  lapply(1:length(data), function(x) {
    lapply(1:length(data[[x]]), function(y) {
      out <- bed[bed[,4] %in% rownames(data[[x]][[y]]),]
      write.table(out, file=paste(fname, names(data[[x]])[y], paste("q", names(data)[x], sep=""), sep="_"),
                  quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
    })
  })  
}

prepForHeatmap <- function(mat, var_ind=c(1:3), N=1000) {
  mat <- apply(mat, 2, pseudoCountNorm)
  mat <- na.omit(mat)
  #mat <- log(mat, 2)
  mat <- sqrt(mat)
  mat_extremes <- quantile(mat, probs=c(.02,.98))
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
