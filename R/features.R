
library("RColorBrewer")

source("~/src/seqAnalysis/R/paths.R")
source("~/src/seqAnalysis/R/blocks.R")
source("~/src/seqAnalysis/R/profiles2.R")
source("~/src/seqAnalysis/R/boot.R")
source("~/src/seqAnalysis/R/modeling.R")



registerDoMC(cores=3)

feature.path <- "~/lib/features_general"
feature_norm_path <- "~/s2/analysis/features/norm"

features_toplot <- factor(1:11, labels=c("refgene_1to3kb_up_chr", "Refgene_1kb_up_chr", "cgi_chr", "Refgene_5_UTR_chr",
                                  "Refgene_CDS_chr", "Refgene_3_UTR_chr", "Refgene_exons_chr",
                                  "Refgene_intron_chr",
                                  "phastCons30way_intergenic_sorted_merge500_thresh500_inter_omp_h3k4me1_default.bed_chr",
                                  "phastCons30way_intergenic_merge500_thresh500_chr",
                                  "intergenic_sub_rmsk_chr"))

features_merge_toplot <- factor(1:11, labels=c("refgene_1to3kb_up_merged",
                                       "Refgene_1kb_up_merged",
                                       "cgi_merged",
                                       "Refgene_5_UTR_merged",
                                       "Refgene_CDS_merged",
                                       "Refgene_3_UTR_merged",
                                       "Refgene_exons_merged",
                                       "Refgene_intron_merged",
                                       "phastCons30way_intergenic_sorted_merge500_thresh500_inter_omp_h3k4me1_default.bed_merged",
                                       "phastCons30way_intergenic_merge500_thresh500_merged",
                                       "intergenic_sub_rmsk_merged"))

features_merge_toplot_short <- factor(1:11, labels=c("1 to 3 kb\n upstream",
                                 "1kb upstream",
                                 "CGI",
                                 "5' UTR",
                                 "CDS",
                                 "3' UTR",
                                 "Exon",
                                 "Intron",
                                 "mOSN\n enhancer",
                                 "Intergenic\n conserved",
                                 "Intergenic"))

features_rmsk <- factor(1:10, labels=c("rmsk_LTR_chr", "rmsk_LINE_chr", "rmsk_SINE_chr",
                               "rmsk_DNA_chr", "rmsk_RNA_chr",
                                "rmsk_rRNA_chr", "rmsk_snRNA_chr", "rmsk_tRNA_chr", "rmsk_srpRNA_chr",
                               "rmsk_Satellite_chr"))
features_rmsk_short <- c("LTR", "LINE", "SINE", "DNA", "RNA", "rRNA", "snRNA", "tRNA","srpRNA", "Satellite")

features_merged_rmsk <- factor(1:10, labels=c("rmsk_LTR_merged", "rmsk_LINE_merged", "rmsk_SINE_merged",
                               "rmsk_DNA_merged", "rmsk_RNA_merged",
                                "rmsk_rRNA_merged", "rmsk_snRNA_merged", "rmsk_tRNA_merged", "rmsk_srpRNA_merged",
                               "rmsk_Satellite_merged"))

features_toplot_nochr <- factor(1:12, labels=c("refgene_1to3kb_up", "cgi", "Refgene_5_UTR",
                                  "Refgene_CDS", "Refgene_3_UTR", 
                                  "Refgene_intron",
                                  "phastCons30way_intergenic_sorted_merge500_thresh500_inter_omp_h3k4me1_default.bed",
                                  "phastCons30way_intergenic_merge500_thresh500",
                                  "rmsk_LINE", "rmsk_SINE", "rmsk_Satellite",
                                  "intergenic_sub_rmsk"))

feature_toplot2 <- factor(1:9, labels=c("refgene_1to3kb_up_chr", "cgi_chr", "Refgene_5_UTR_chr",
                                  "Refgene_CDS_chr", "Refgene_3_UTR_chr", "Refgene_exons_chr",
                                  "Refgene_intron_chr", "omp_mk4_intergenic_inter_cons.bed_chr",
                                  "phastCons30way_intergenic_merge500_thresh500_chr"))

featureScaleFill5 <- scale_fill_manual(values=brewer.pal(5, "Set1"))
featureScaleFill4 <- scale_fill_manual(values=brewer.pal(4, "Set1"))
featureScaleFill3 <- scale_fill_manual(values=brewer.pal(3, "Set1"))


featureClasses <- factor(1:4, labels=c("Upstream", "Exons", "Transcript", "Intergenic"))
featureMatrix <- data.frame(features_toplot, features_merge_toplot, features_merge_toplot_short, featureClasses[c(1,1,1,2,2,2,3,3,4,4,4)])
featureClasses5 <- factor(1:5, labels=c("Upstream", "Exons", "Transcript", "Intergenic", "Repeat"))

enhancerClasses <- data.frame(cluster=1:14, class=factor(1:4, labels=c("5hmC variable", "5hmC/5mC variable", "5mC variable", "Invariant"))[c(1,1,1,1,2,2,2,2,3,4,3,3,3,3)])

## Make summaries --------------------------
makeFeatureMatrix2 <- function(feature, set = "all", select=NULL, data_type, transf=NULL, write=TRUE) {
  if (set=="cells_norm" | set=="cells") {
    samples <- samples.cells_norm
  } else if (set=="unnorm") {
    samples <- samples.cells_raw
  } else if (set=="cells_rpkm") {
    samples <- samples.cells_rpkm
  } else if (set=="cells_rpkm2") {
    samples <- c("omp_hmc_120424_rpkm", "ngn_hmc_rpkm", "icam_hmc_rpkm",
                 "omp_mc_rpkm", "ngn_mc_rpkm", "icam_mc_rpkm")
  } else if (set=="cells_strand") {
    samples <- c("omp_hmc_120424_rmdup_plus_sub_minus", "ngn_hmc_120424_rmdup_plus_sub_minus",
                 "icam_hmc_120424_rmdup_plus_sub_minus", "omp_mc_rmdup_plus_sub_minus",
                 "ngn_mc_rmdup_plus_sub_minus", "icam_mc_rmdup_plus_sub_minus")
  } else if (set=="cells_noM") {
    samples <- samples.cells_noM
  } else if (set=="cells_full") {
    samples <- c("omp_hmc_120424_full", "ngn_hmc_120424_full", "icam_hmc_120424_full",
                 "omp_mc_full", "ngn_mc_full", "icam_mc_full")
  } else if (set=="cells_bam") {
    samples <- c("omp_hmc_120424_rmdup", "ngn_hmc_120424_rmdup", "icam_hmc_120424_rmdup")
  } else if (set=="medips_rf_1" | set=="medips_rf_2") {
    samples <- samples.medips_rf
  } else if (set=="rpm_avg_2") {
    samples <- samples.cells_norm
  } else if (set=="tfo_unnorm") {
    #data_type <- "unnorm_2"
    samples <- samples.tfo_unnorm
  } else if (set=="tfo") {
    samples <- samples.tfo
  } else if (set=="d3a") {
    samples <- samples.d3a_norm
  } else if (set=="d3a_unnorm") {
    samples <- samples.d3a_norm
  } else if (set=="d3a_rpkm") {
    samples <- samples.d3a_rpkm
  } else if (set=="d3a_2") {
    samples <- samples.d3a_2
  } else if (set=="d3a_strand") {
    samples <- c("moe_d3a_wt_hmc_plus_sub_minus", "moe_d3a_ko_hmc_plus_sub_minus",
                 "moe_d3a_wt_mc_plus_sub_minus", "moe_d3a_ko_mc_plus_sub_minus",
                 "moe_d3a_wt_in_plus_sub_minus")
  } else if (set=="all") {
    samples <- list.files(paste(feature_norm_path, feature, sep="/")) 
  } else if (set=="d3a_mrna") {
    samples <- c("moe_d3a_wt_mrna_rpkm", "moe_d3a_ko_mrna_rpkm")
  } else if (set=="d3a_mrna_2") {
    samples <- c("moe_d3a_wt_mrna", "moe_d3a_ko_mrna")
  } else if (set=="cells_nuc") {
    samples <- c("omp_nuc_0123", "icam_nuc_01234")
  } else if (set=="d3a_nuc") {
    samples <- samples.d3a_nuc
  } else if (set=="d3a_nuc_sub") {
    samples <- samples.d3a_nuc_sub
  } else if (set=="d3a_kde") {
    samples <- c("d3xog_wt_nuc_478_p1", "d3xog_ko_nuc_256_p1")
  } else if (set=="tt3_min") {
    #samples <- c("omp_hmc_120424_rpkm", "ott3_1_hmc_rpkm", "ott3_2_hmc_rpkm", "ott3_2_hmc_21M_rpkm", "ngn_hmc_rpkm", "icam_hmc_rpkm",
    #             "omp_mc_rpkm", "ott3_1_mc_rpkm", "ott3_2_mc_rpkm", "ott3_2_mc_rmdup_17M_rpkm", "ngn_mc_rpkm", "icam_mc_rpkm")
    samples <- c("omp_hmc_120424_rpkm", "ott3_1_hmc_rpkm", "ott3_2_hmc_rpkm",
                 "omp_mc_rpkm", "ott3_1_mc_rpkm", "ott3_2_mc_rpkm")
  } else if (set=="tt3_sub") {
    samples <- c("ott3_2_hmc_21M_rpkm_sub_omp_hmc_120424_rpkm")
  } else if (set=="cpg") {
    samples <- c("mm9")
  } else if (set=="encode_dnase") {
    samples <- c("wgEncodeUwDnaseCerebrumC57bl6MAdult8wksAlnRep1", "wgEncodeUwDnaseCerebellumC57bl6MAdult8wksAlnRep1",
                 "wgEncodeUwDnaseHeartC57bl6MAdult8wksAlnRep1", "wgEncodeUwDnaseLiverC57bl6MAdult8wksAlnRep1",
                 "wgEncodeUwDnaseRetinaC57bl6MAdult1wksAlnRep1")
  } else if (set=="encode_mk4") {
      samples <- c("encode_cbellum_h3k4me1", "encode_cortex_h3k4me1", "encode_heart_h3k4me1",
                   "encode_liver_h3k4me1", "encode_lung_h3k4me1")
  }
  if (!is.null(select)) {
    samples <- sapply(samples, function(x) paste(x, select, sep="_"))
  }
  print(samples)
  vals <- foreach (sample=samples, .combine="cbind") %dopar% {
    print(paste(feature_norm_path, data_type, feature, sample, sep="/"))
    data <- read.delim(paste(feature_norm_path, data_type, feature, sample, sep="/"), header=FALSE)
    val <- data[,ncol(data)]
    if (!is.null(transf)) val <- do.call(transf, list(val))
    names(val) <- data[,4]
    return(val)
  }
  vals <- as.matrix(vals)
  colnames(vals) <- samples
  out_path <- paste(feature_norm_path, data_type, "summaries", sep="/")
  if (!file.exists(out_path)) dir.create(out_path)
  fname <- paste(feature_norm_path, data_type, "summaries", paste(set, feature, sep="_"), sep="/")
  if (!is.null(transf)) fname <- paste(fname, transf, sep="_")
  if (write) write.table(vals, file=fname, quote=FALSE, sep="\t")
  return(vals)
}

makeFeatureMatrix2.all <- function(set, data_type, transf=NULL) {
  files <- list.files(feature.path)
  files <- files[grep("chr", files)]
  for(file in files) {
    print(file)
    if (file.exists(paste(feature_norm_path, data_type, "summaries", paste(set, file, transf, sep="_"), sep="/"))) {
      print("File exists")
      next
    }
    tryCatch(makeFeatureMatrix2(file, set=set, data_type=data_type, transf=transf, write=TRUE), error = function(e) {
            print(paste("Skipping", file, sep=" "))
                  print(e)
                  return
          })
  }
}


# Summary functions -------------------------------------------------------


statCollect <- function(set, data_type, feature, transf_name, transf=NULL) {
  data_path <- paste(feature_norm_path, data_type, "summaries",
                           paste(set, feature, sep="_"), sep="/")
  if (!is.null(transf_name)) data_path <- paste(data_path, transf_name, sep="_")
  print(data_path)
  data <- read.delim(data_path, header=TRUE, row.names=NULL)
  data_melt <- melt(data)
  colnames(data_melt)[1] <- "obs"
  data_melt$feature <- feature
  labels <- sapply(colnames(data)[2:ncol(data)], str_split, "_")
  celltype <- unique(unlist(lapply(labels, function(x) x[1])))
  ip <- unique(unlist(lapply(labels, function(x) x[2])))
  data_melt$celltype <- rep(celltype, each=nrow(data))
  data_melt$ip <- rep(ip, each=nrow(data) * length(celltype))
  if (!is.null(transf)) data_melt <- transform(data_melt, value=do.call(transf, list(value)))
  return(data_melt)
}

statSummary <- function(set, data_type, feature, transf_name=NULL, transf=NULL, FUN, ...) {
  data_path <- paste(feature_norm_path, data_type, "summaries",
                           paste(set, feature, sep="_"), sep="/")
  if (!is.null(transf_name)) data_path <- paste(data_path, transf_name, sep="_")
  data <- read.delim(data_path, header=TRUE, row.names=NULL)
  #return(data)
  if (!is.null(transf)) {
    if (ncol(data) == 2) {
      data[,2] <- do.call(transf, list(data[,2]))
    } else {
      data[,2:ncol(data)] <- apply(data[, 2:ncol(data)], 2, function(x) do.call(transf, list(x)))
    }  
  }
  if (ncol(data) == 2) {
    return(do.call(FUN, list(data[,2])))
  } else { 
    return(apply(data[,2:ncol(data)], 2, FUN, ...))
  }  
}

statSummary.all <- function(set, data_type, toplot="general", action="summary", transf_name="sqrt", transf=NULL, FUN=mean, ...) {
    features <- ""
    if (toplot=="general") {
      features <- as.character(features_toplot)
    } else if (toplot=="rmsk") {
      features <- as.character(features_rmsk)
    } else {

  
      features <- list.files(paste(feature_norm_path, data_type, "summaries", sep="/"))
                                        #print(features)
      features <- lapply(features, str_split, paste(set, "_", sep=""))
                                        #print(features)
      features <- na.omit(unlist(lapply(features, function(x) x[[1]][2])))
                                        #print(features)
    }
  
  out <- foreach (feature=features, .combine="rbind") %dopar% {
    print(feature)
    if (action=="summary") {
      result <- statSummary(set, data_type, feature, transf_name=transf_name, transf=transf, FUN, ...)
    } else if (action=="collect") {
      result <- statCollect(set, data_type, feature, transf_name=transf_name)
    }
    return(result)
  }
  if (action=="summary") rownames(out) <- features
  return(out)
}

statSummary.test <- function(summary, type="wilcox", col_split="feature", col_compare="celltype", subsample=0) {
  registerDoMC(cores=2)
  summary_split <- split(summary, summary[,col_split])
  samples <- as.character(unique(summary[, col_compare]))
  print(samples)
  samples_comb <- combn(samples, 2)
  samples_comb_name <- unlist(apply(samples_comb, 2, paste, collapse="_"))
  pvalues <- foreach(group=summary_split) %dopar% {

    pvalue <- apply(samples_comb, 2, function(comb) {
      if (type=="wilcox") {
        val <- 0
        if (subsample == 0) {
          val <- wilcox.test(group$value[group[,col_compare]==comb[1]],
                                group$value[group[,col_compare]==comb[2]])$p.value
          return(val)
        } else {
          val <- sapply(1:100, function(x) {
            obs <- sample(unique(group$obs), subsample)
            group_sub <- group[group$obs %in% obs,]
            wilcox.test(group_sub$value[group_sub[,col_compare]==comb[1]],
                        group_sub$value[group_sub[,col_compare]==comb[2]])$p.value
          })
          return(mean(val))
        }
      } else if (type=="perm") {
        with(feature, permutationTest(value[celltype==comb[1]], value[celltype==comb[2]]))
      }  
    })
    pvalue <- unlist(pvalue)
    names(pvalue) <- samples_comb_name
    a <- gc()
    return(pvalue)
  }
  names(pvalues) <- names(summary_split)
  return(pvalues) 
}
  
statSummary.allCI <- function(set, data_type, transf=NULL, boot=TRUE, ...) {
  stat <- melt(statSummary.all(set=set, data_type=data_type, transf=transf, FUN=mean, ...))
  if (boot) {
    CI <- bootCI
  } else {
    CI <- ci95
  }  
  CIup <- melt(statSummary.all(set=set, data_type=data_type, transf=transf,
                               FUN=CI, stat=boot.samplemean, bound="upper" ))
  CIlow <- melt(statSummary.all(set=set, data_type=data_type, transf=transf,
                                FUN=CI, stat=boot.samplemean, bound="lower"))
  colnames(stat) <- c("feature", "sample", "median")
  colnames(CIup) <- c("feature", "sample", "CIup")
  colnames(CIlow) <- c("feature", "sample", "CIlow")
  stat_ci <- merge(stat, CIup, all.y=TRUE)
  stat_ci <- merge(stat_ci, CIlow, all.y=TRUE)
  return(stat_ci)
}

statSummary.allSEM <- function(set, data_type, FUN=sem) {
  stat <- melt(statSummary.all(set=set, data_type=data_type, FUN=median))
  SEM <- melt(statSummary.all(set=set, data_type=data_type, FUN=FUN))
  colnames(stat) <- c("feature", "sample", "median")
  colnames(SEM) <- c("feature", "sample", "SEM")
  return(merge(stat, SEM, all.y=TRUE))
}

statSummary.allSD <- function(set, data_type, trim=.05, ...) {
  stat <- melt(statSummary.all(set=set, data_type=data_type, FUN=mean, trim=trim))
  sd <- melt(statSummary.all(set=set, data_type=data_type, FUN=sd_trim, trim=trim))
  colnames(stat) <- c("feature", "sample", "median")
  colnames(sd) <- c("feature", "sample", "sd")
  return(merge(stat, sd, all.y=TRUE))
}

statSummary.allIQR <- function(set, data_type, toplot="general", transf_name="sqrt", transf=NULL, ...) {
  stat <- melt(statSummary.all(set=set, data_type=data_type, toplot=toplot, transf_name=transf_name, transf=transf, FUN=median, ...))
 
  up <- melt(statSummary.all(set=set, data_type=data_type, toplot=toplot, transf_name=transf_name, transf=transf,
                               FUN=iqr, bound="upper" ))
  low <- melt(statSummary.all(set=set, data_type=data_type, toplot=toplot,
                              transf_name=transf_name, transf=transf,
                                FUN=iqr, bound="lower"))
  colnames(stat) <- c("feature", "sample", "median")
  colnames(up) <- c("feature", "sample", "up")
  colnames(low) <- c("feature", "sample", "low")
  stat_ci <- merge(stat, up, all.y=TRUE)
  stat_ci <- merge(stat_ci, low, all.y=TRUE)
  return(stat_ci)
}

## Normalize values of feature matrix by median value of specified feature
statSummary.allNorm <- function(set, data_type, toplot=TRUE, transf_name=NULL, transf=NULL, norm_col="feature", norm="intergenic_sub_rmsk_chr") {
  data <- statSummary.all(set=set, data_type=data_type, toplot=toplot, action="collect", transf_name=transf_name, transf=transf)
#  data.norm <- data[data$feature == norm,]
    data.norm <- data[data[,norm_col] == norm,]
  data.norm.median <- ddply(data.norm, c("variable", "ip", norm_col), summarize, value.median=median(value))
  norm_feature_name <- paste(norm, "median", sep="_")
  data.norm.median <- data.frame(obs="X",
                                 variable=data.norm.median$variable,
                                 value=data.norm.median$value.median,
                                 feature=norm_feature_name,
                                 celltype=data.norm.median$celltype,
                                 ip=data.norm.median$ip)
  data <- rbind(data, data.norm.median)
  print("start norm")
  data.split <- split(data, data$variable)
  data.op <- lapply(data.split, function(x) {x$value = x$value/x$value[x$feature==norm_feature_name]
                                           return(x)})
  data <- do.call("rbind", data.op)
  data <- ddply(data, .(variable, feature, celltype, ip), summarize, value.median=median(value), up=iqr(value, bound="upper"), low=iqr(value, bound="lower"))
  return(data)
}



# Intersection functions --------------------------------------------------


processIntersectSummary <- function(summary, features=features_merge_toplot) {
  data <- read.delim(summary)
  data$internal_norm <- with(data, fraction_norm/sum(fraction_norm))
  data$feature.factor <- features[match(data$feature, as.character(features))]
  data$feature.pretty <- features_merge_toplot_short[as.numeric(data$feature.factor)]
  data$class <- featureMatrix[match(data$feature.pretty, featureMatrix[,3]),4]
  return(data)
}


processIntersectSummary.batch <- function(set, features=features_merge_toplot) {
  path <- "~/s2/data/homer/peaks/intersections"
  set_path <- paste(path, set, sep="/")
  peak_sets <- list.files(set_path)
  out <- foreach(peak_set=peak_sets, .combine="rbind") %do% {
    summary_file <- paste(set_path, peak_set, "summary", sep="/")
    data <- processIntersectSummary(summary_file, features)
    data$peak_set <- peak_set
    return(data)
  }
  return(out)
}


# Utility functions -------------------------------------------------------


sem <- function(vals) {
  return(sd(vals)/(length(vals)))
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

feature_wilcox <- function(test, ref) {
  return(wilcox.test(test, ref)$p.value)
}

feature_ks_batch <- function(collect_df, key_col, reference_key) {
  sub <- ddply(collect_df, c("feature", "ip"), function(values, ind) {cell_compare(values, ind, key_col, reference_key)})
  return(sub)
}


#! not functional
cell_compare <- function(values, ind, key_col, reference_key) {
  print(dim(values))
  #print(names(values))
  #sub_values <- values[,ind]
  #print(dim(sub_values))
  keys <- unique(values[,key_col])
  sub <- unlist(lapply(keys, function(key) ks_pvalue(values[values[,key_col]==key,"value"],
                                                     values[values[,key_col]==reference_key, "value"])))
  #sub <- unlist(lapply(keys, function(key) permutationTest(values[values[,key_col]==key,"value"],
  #                                                   values[values[,key_col]==reference_key, "value"], FUN="median")))
  names(sub) <- keys
  #sub <- ddply(values, c("celltype"), summarize,
  #             ks_pvalue=ks_pvalue(value, values[values$key_col==reference_key,value]))
  return(sub)
}

ks_pvalue <- function(data1, data2) {
  return(wilcox.test(data1, data2)$p.value)
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

splitByQ <- function(data, column, q=c(0, .25, .5, .75, 1), fname=NULL) {
  qs <- quantile(data[,column], q)
  data_cut <- cut(data[,column], breaks=qs, labels=FALSE)
  data_out <- data.frame(rownames(data), data_cut)
  if (!is.null(fname)) write.table(data_out, file=fname, quote=FALSE, sep="\t", row.names=F, col.names=F)
  return(data_out)
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

makeFeatureDF <- function(set="cells", data_type="raw", FUN=mean) {
  files <- list.files(paste(feature2.path, "summaries", sep="/"))
  files <- files[grep(set, files)]
  files <- files[grep(data_type, files)]
  files_names <- lapply(files, str_split, "_")
  print(files)
  #return(files_names)
  files_names <- lapply(files_names, function(x) {
    #print(x)
    sel <- x[[1]][-grep(set, x[[1]])]
    #return(sel)
    sel <- sel[-grep(data_type, sel)]
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
  stat.df$data_type <- "stat"
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

computeScoreRatios <- function(v1, v2) {
  return(log((v1+1)/(v2+1))*v1)
}

compareFeatures <- function(feature, sample_list, value_type="raw") {
  data <- lapply(sample_list, function(sample) {
    read.delim(paste(feature2.path, sample, feature, sep="/"), header=FALSE)
  })
  
}

ecdf <- function(d) {
  h <- hist(d, breaks=50, plot=FALSE)
  counts <- h$counts
  count_sum <- sum(counts)
  vals <- vector('numeric', length=length(counts))
  vals[1] <- counts[1]
  for (i in 2:length(vals)) {
     vals[i] <- counts[i] + vals[i-1]
  }
  vals <- vals / count_sum
  return(cbind(breaks=h$mids, vals=vals))
}

