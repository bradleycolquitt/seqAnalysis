library(stats)
library(itertools)
library(foreach)
library(multicore)
library(mclust)

#source("~/src/R/LEA/dipdata.R")
source("~/src/R/LEA/seqdata.R")
source("~/src/R/plotUtil.R")
source("~/src/seqAnalysis/R/profiles2.R")
source("~/src/seqAnalysis/R/util.R")

registerDoMC(cores=6)
corrs <- function(data1, data2, method="spearman", N) {
	index <- seq(1, length(data1) - N , by = N)
	cors <- c()
	for (i in index){
		sub1 <- data1[seq(i, i + N)]
                sub2 <- data2[seq(i, i + N)]
		cors <- c(cors,cor(sub1, sub2, method=method))
		}
	return(cors)
	}
	
corr.null<-function(data1, data2, n){
	index<-seq(1,nrow(data)-n,by=n)
	cors.n<-c()
	for (i in index){
		r<-seq(i,i+n)
		r.rand1<-sample(r,n)
		r.rand2<-sample(r,n)
		cors.n<-c(cors.n,cor(data1[r.rand1],data2[r.rand2]))
		}
	cors.n
	}

corrMedips <- function(data1, data2, data_type="norm") {
  read_vals.1 <- getDataByType(data1, data_type)
  read_vals.2 <- getDataByType(data2, data_type)

  sample_corr <- corrs(read_vals.1, read_vals.2, 100)
  sapmle_cors_null <- corrs.null(read_vals.1, read_vals.2, 100)

  ks <- ks.text(sample_corr, sample_corr_null)
  print(ks)
  return(mean(sample_corr, na.rm=T))
}

corrDipData <- function(data1, data2, data_type="norm", method="spearman",
                        N=1000, test=FALSE) {
  dd1 <- load.DipData(data1)
  dd2 <- load.DipData(data2)
  
  dd1.data <- dd1$genome.data[dd1$genome.data[,1]!=22, data_type]
  dd2.data <- dd2$genome.data[dd1$genome.data[,1]!=22, data_type]
  sample_corr <- corrs(dd1.data, dd2.data, method=method, N=N)
  if(test) {
    sample_corr_null <- corrs.null(dd1.data, dd2.data, N)
    ks <- ks.test(sample_corr, sample_corr_null)
    print(ks)
  }
  
  return(mean(sample_corr, na.rm=TRUE))
}

meanWindows <- function(data, step) {
  end <- length(data)/step
  ind <- 0
  mod <- length(data)%%step
  if (mod == 0) {
    ind <- rep(1:end, each=step)
  } else {
    ind <- c(rep(1:end, each=step), rep(end+1, times=mod)) 
  }
  return(tapply(data, ind, mean))
}
  
corrWig <- function(wig1, wig2, method="spearman", N=1000, step=1) {
  
  wig1_chr <- list.files(wig1)
  wig2_chr <- list.files(wig2)
  #print(wig1_chr)
  common_chr <- wig1_chr[wig1_chr %in% wig2_chr]
  filter <- c(grep("random", common_chr), grep("NT", common_chr))
  if (length(filter) > 0) common_chr <- common_chr[-filter]
  corrs <- foreach(chr=common_chr, .combine="c", .verbose=FALSE) %dopar% {
    #print(chr)
    data1 <- scan(paste(wig1, chr, sep="/"), skip=1, quiet=TRUE)
    data2 <- scan(paste(wig2, chr, sep="/"), skip=1, quiet=TRUE)
    if (step > 1) {
      #print(paste("Binning ", chr, sep=""))
      time <- Sys.time()
      data1 <- meanWindows(data1, step)
      data2 <- meanWindows(data2, step)
      #print(paste(chr, Sys.time() - time, sep=" "))
    }
    if (length(data1) != length(data2)) {
      stop(paste("Non-matching wig lengths: ", chr, sep=""))
    }
    tryCatch(cor(data1, data2, method=method), error=function(e) return(NA))   
  }
  return(list(mean=mean(corrs), sd=sd(corrs), se=sd(corrs)/length(corrs)))
}

stepCorrWig <- function(wig1, wig2, method="spearman", steps) {
  step_corrs <- list()
  for (step in steps) {
    print(step)
    step_corrs <- c(step_corrs, list(corrWig(wig1, wig2, method=method, step=step)))
  }
  names(step_corrs) <- steps
  return(step_corrs)
}

hmedips <- paste(c("omp", "ngn","icam"), "hmedip", sep="_")
medips <- paste(c("omp", "ngn","icam"), "medip", sep="_")
mk4 <- paste(c("omp", "ngn"), "mk4", sep="_")
corrDNAmodAndmk4 <- function(data_type="raw", method="spearman", N=100, fname=NULL) {
  #print(mk4)
  dipdata <- foreach(data=c(hmedips, medips), .combine="cbind") %do% {
    
    dd <- load.DipData(data)    
    genome.data <- dd$genome.data
    dd.data <- genome.data[!genome.data[,1]==22, data_type]
    if (data_type == "raw") dd.data <- (dd.data * 1e6) / sum(dd.data)
    return(dd.data) 
  }
  
  seqdata <- foreach(data=mk4, .combine="cbind") %do% {
    cat("mk4\n")
    print(data)
    sd <- load.SeqData(data)
    print("loaded")
    genome.data <- sd$genome.data
    sd.data <- genome.data[!genome.data[,1]==22, "rpm"]   
    return(sd.data) 
  }
  
  mat <- cbind(dipdata, seqdata)
  #return(mat)
  colnames(mat) <- c(hmedips, medips, mk4)
  #return(mat)
  corr_mat <- corrMatrix(mat, method=method, N=N, fname=fname)
  return(corr_mat)
}



#Combinatorial patterns of histone acetylations and methylations in the human genome
corrMatrix <- function(set, method, N, fname) { 
 # cat("Getting data...\n")
 # set_data <- foreach(data=set, .combine="cbind") %do% {
 #   dd <- load.DipData(data)    
 #   genome.data <- dd$genome.data
 #   dd.data <- genome.data[!genome.data[,1]==22, data_type]   
 #   return(dd.data) 
 # }
  #return(set)
 
  cat("Removing zeros...\n")
  set_names <- colnames(set)
  set <- removeZeros(set)
  gc()
  cat("Correlating...\n")
  comb <- combn(set_names, 2)
  registerDoMC(cores=6)
  corr_vals <- foreach(pair=iter(comb, by="col"), .inorder=TRUE, .combine="c") %dopar% {
    cat(paste(pair[1], " and ", pair[2], "\n", sep=""))
    val <- corrs(set[,pair[1]], set[,pair[2]], method=method, N=N)
    return(mean(val, na.rm=TRUE))
  }
  
  mat <- fillMatrix(set_names, corr_vals)
  dimnames(mat) <- list(set_names, set_names)
  if (!is.null(fname)) write.table(mat, file=paste("~/storage/analysis/cor/", fname, sep=""),
              quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
  return(mat)
}

mclustDNAmodRNA <- function(feature, dna=NULL, set, value_type="", rna=NULL) {
  #rna.path <- "~/data/rna/omp_ngn_icam_rna.txt"
  #rna.path <- "~/data/rna/cv_iv_cd_id_rna"
  #rna.path <- "~/storage/data/rna/cuffdiff/moe_wt_mrna_moe_d3a_mrna/gene_exp.diff"
  #rna.path <- "~/s2/data/rna/cuffdiff/omp_mrna_ngn_mrna/gene_exp.diff"
  rna.path <- "~/s2/analysis/rna/summaries/omp_ngn_icam_mrna_dup_nozero_log2"
  #dna.path <- paste("~/analysis/mprofiles/features/feature_summaries/", paste(dip, "refgene_long", sep="_"), sep="")
  dna.path <- paste("~/s2/analysis/features/norm", value_type, "summaries", paste(set, feature, sep="_"), sep="/")
  #dna.path <- "~/storage/analysis/mprofiles/features/moe_wt_hmc_moe_dnmt3a_hmc/refgene_noclust"
  rna.data.full <- read.delim(rna.path)
  #rna.data <- matrix(c(rna.data.full$value_1, rna.data.full$value_2), nrow=nrow(rna.data.full), ncol=2,

  #rna.data <- matrix(c(rna.data.full$value_1, rna.data.full$value_2), nrow=nrow(rna.data.full), ncol=2,
  #                   dimnames=list(rna.data.full$test_id, c("mOSN", "GBC")))
                                        #dimnames(rna.data) <- list(rna.data.full$test_id, c("WT", "Dnmt3a"))
  #rna.data <- removeZeros(rna.data)
  dna.data <- read.delim(dna.path)
  dna.data <- removeZeros(dna.data)
  #dna.data <- apply(dna.data, 2, pseudoCountNorm)
  #return(dna.data)
                                        #dimnames(dna.data) <- list(rownames(dna.data.full), c("WT", "Dnmt3a"))
  #rna.sample.data <- log(rna.data[, grep(rna, colnames(rna.data))],2)
  rna.sample.data <- rna.data.full[,grep(rna, colnames(rna.data.full))]
  #dna.sample.data <- log(dna.data[, grep(dna, colnames(dna.data))], 2)
  dna.sample.data <- sqrt(dna.data[, grep(dna, colnames(dna.data))])
  rna.select.data <- rna.sample.data[match(rownames(dna.data), rownames(rna.data.full))]
  #rna.select.data <- rna.data[match(rownames(dna.data), rownames(rna.data))]
  sample.data <- na.omit(cbind(rna.select.data, dna.sample.data))
#  return(sample.data)
  sample.mc <- Mclust(sample.data, G=2)
  return(list(sample.mc, sample.data))
  #corr <- cor(dna.sample.data, rna.select.data, method="spearman", use="complete.obs") 
  #col <- rgb(0,0,0,.2)
  #x11()
  #plot(sample.data[,2], sample.data[,1], col=col, pch=20, ylim=c(400,1000))
  #return(corr)
}

mclust.several <- function(feature, set="cells", value_type="") {
  if (set=="cells") {
    dna.samples <- samples.cells_norm
    rna.samples <- c("omp", "ngn", "icam")
  } else if (set=="medips_rf_1" | set=="medips_rf_2") {
    dna.samples <- samples.medips_rf
    
  } else if (set=="unnorm") {
    dna.samples <- samples.cells_norm
    rna.samples <- c("omp", "ngn", "icam")
  } else if (set=="unnorm_2") {
    dna.samples <- samples.cells_norm
    rna.samples <- c("omp", "ngn", "icam")
  }

  samples <- list()
  samples_paste <- list()
  for (i in 1:length(dna.samples)) {
    for(j in 1:length(rna.samples)) {
      samples <- c(samples, list(c(dna.samples[i], rna.samples[j])))
      samples_paste <- c(samples_paste, list(paste(c(dna.samples[i], rna.samples[j]), collapse="_")))
    }             
  }
  #return(samples)
  registerDoMC(cores=6)
  mcs <- foreach(sample=samples, .inorder=TRUE) %dopar% {
    print(sample)
    mc <- mclustDNAmodRNA(feature=feature, dna=sample[1], set=set, value_type=value_type, rna=sample[2])
    gc()
    return(mc)
  }
  names(mcs) <- unlist(samples_paste)
  return(mcs)
}
mclustSplitClass <- function(mclust, data, class=2, plot=FALSE, corr=TRUE) {
  classif <- mclust$classification
  data.class <- data[classif==class,]
  if (plot) plot(data.class[,1], data.class[,2])
  if (!corr) return(data.class)
  data.cor <- cor(data.class[,1], data.class[,2], method="spearman")
                                        #data.cor <- cor.test(data.class[,1], data.class[,2], method="spearman")
  #data.cor <- lm(data.class[,2]~ data.class[,1])
  #data.cor <- glm(data.class[,2] ~ data.class[,1])
  return(data.cor)
}

mclustFit <- function(mclust, data, class=2) {
  data.class <- NULL
  if (class > 0) {data.class <- mclustSplitClass(mclust, data, class)
  } else {data.class <- data[data[,1] >= -1,]}
  data.fit <- lm(data.class[,2]~data.class[,1])
  return(data.fit)
}

cust.topo <- topo.colors(50)
cust.topo[1] <- "white"

mclustPlotSurface <- function(mclust, data, type="image", x.lim=c(-8, 8), y.lim=c(0, 1), fname=NULL, corr.val=NULL, corr.val.pos=c(0,0)) {
  # if (is.null(fname)) {x11()
  # } else {
  #   pdf(file=paste("~/s2/analysis/features/plots",
  #         paste(fname, ".pdf", sep=""), sep="/"),
  #       width=5, height=5)
  # }                     
   #y.lim <- c(min(data[,2]), 0)
   #y.lim <- c(-6, -2)           
   surfacePlot(data=data, what="density", type=type,
               parameters=mclust$parameters,
               col=topo.colors(100),
               ylim=y.lim, xlim=x.lim,
               cex.axis=1.2,
               ann=FALSE)
   #mtext("log2 FPKM", side=1, line=2.5, cex=1.2)
   #mtext(expression(sqrt("RPM")), side=2, line=2.5, cex=1.2)
   if (!is.null(corr.val)) {
     corr.val <- round(corr.val, digits=2)
     text(corr.val.pos[1], corr.val.pos[2], paste(expression(rho), " = ", as.character(corr.val), sep=""), col="white", pos=1, cex=2)
   }
   if(!is.null(fname)) dev.off()
 }

mclustComplete <- function(dna, dip, rna) {
  data <- plotDNAmodRNA(dna, dip, rna)
  mclustPlotSurface(data[[1]], data[[2]])
}

sectionByFullQ <- function(rna_file, dna_file, Q) {
  rna <- read.delim("~/data/rna/omp_ngn_icam_rna.txt")
  dna <- read.delim(paste("~/storage/analysis/mprofiles/features/feature_summaries/", dna_file, sep=""))
  dna.unwrap <- c(dna[,1], dna[,2], dna[,3])
  q.val <- quantile(dna.unwrap, probs=Q)
  dna.filter <- list()
  
  for (i in 1:ncol(dna)){
    ind <- dna[,i] >= q.val
    m <- match(rownames(dna[ind,]), rownames(rna))
    vals <- cbind(dna[ind,i], rna[m,i])
    #dna.filter[[i]] <- dna[ind, i]
    #names(dna.filter[[i]]) <- rownames(dna[ind,])
    
    #dna.filter[[i]] <- cbind(dna.filter[[i]], rna[m, i])
    rownames(vals) <- rownames(dna[ind,])
    colnames(vals) <- c("dna", "rna")
    vals <- na.omit(vals)
    dna.filter[[i]] <- vals
  }
  names(dna.filter) <- colnames(dna)
  return(dna.filter)
}
TopHmcRnaDensity <- function(val_a, val_b, fname=NULL) {
  if (is.null(fname)) {
    X11()
  } else {
    pdf(fname, 5, 5)
  }
  val_a_density <- density(val_a)
  val_b_density <- density(val_b)
  max_y <- max(val_a_density$y)
  max_y <- max(c(max_y, val_b_density$y))
  print(max_y)
  plot(1,1, type="n", xlim=c(-10,10), ylim=c(0, .25), xlab="log2 FPKM", ylab="Density", axes=FALSE)
  lines(val_a_density$x, val_a_density$y, col=col3[3])
  lines(val_b_density$x, val_b_density$y, col=col3[2])
  axis(2, seq(0, .25, .05))
  axis(1)
  #box()
  legend(-8, max_y-(max_y*.1), c("OSN", "GBC"), col=col3[3:2], bty="n", lty=1)
  if (!is.null(fname)) dev.off()
}


fillMatrix <- function(set, vals) {
   mat <- matrix(1, nrow=length(set), ncol=length(set))
   .filler <- function(mat) {
     prev <- 0
     for(i in 1:(length(set) - 1)) {
       step <- prev + length(set) - i
       mat[i, (i+1):ncol(mat)] <- vals[(prev + 1):step]
       prev <- step
    }
    return(mat)
   }
   mat <- .filler(mat)
   mat <- .filler(t(mat))
   return(mat)
 }

distanceMatrix <- function(mat, method="euclidean") {
  mat.nz <- mat[apply(mat, 1, prod) > 0,]
  dist.mat <- dist(t(mat.nz), method=method)
  return(dist.mat)
}

