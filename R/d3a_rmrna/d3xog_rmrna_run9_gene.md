D3xOG rmRNA -- run9
========================================================

Cuffdiff parameters:
* version 2.1.1
* -b /seq/lib/indexes/mm9.fa 
* -M /home/user/lib/rmsk/rmsk.gtf 
* --compatible-hits-norm
* -u
* --library-norm-method classic-fpkm
* --dispersion-method pooled
* genes_chr_protein.gtf

Samples
  * wt
    * 121126/d3xog_wt_rmrna_blank/d3xog_wt_rmrna_blank.bam
  * het 
    * 130326/d3xog_het_rmrna_blank/d3xog_het_rmrna_blank.bam
  * ko
    * 121126/d3xog_ko_rmrna_blank/d3xog_ko_rmrna_blank.bam
    * 130326/d3xog_ko_rmrna_rep2_blank/d3xog_ko_rmrna_rep2_blank.bam
    

```r
opts_chunk$set(warning = FALSE, message = FALSE, error = FALSE)
library(plyr)
```

```
## Attaching package: 'plyr'
```

```
## The following object(s) are masked from '.env':
## 
## unrowname
```

```r
library(reshape2)
library(ggplot2)
library(gridExtra)
```

```
## Loading required package: grid
```

```r
suppressPackageStartupMessages(source("~/src/seqAnalysis/R/features.R"))
suppressPackageStartupMessages(source("~/src/seqAnalysis/R/image.R"))
suppressPackageStartupMessages(source("~/src/seqAnalysis/R/ggplot2.R"))
```


  

```r
samples <- c("omp_hmc_rep1_q30_rmdup_extend300_mean_omp_hmc_rep2_q30_rmdup", 
    "d3xog_het_hmc_sort_q30_rmdup", "d3xog_ko_hmc_sort_q30_rmdup", "omp_mc_rep1_q30_rmdup_extend300", 
    "d3xog_het_mc_sort_q30_rmdup", "d3xog_ko_mc_sort_q30_rmdup", "ngn_hmc_rep1_q30_rmdup_extend300_mean_ngn_hmc_rep2_q30_rmdup", 
    "icam_hmc_rep1_q30_rmdup_extend300_mean_icam_hmc_rep2_q30_rmdup", "ngn_mc_rep1_q30_rmdup_extend300_mean_ngn_mc_rep2_q30_rmdup", 
    "icam_mc_rep1_q30_rmdup_extend300_mean_icam_mc_rep2_q30_rmdup")
data <- lapply(samples, function(x) makeImage(x, "gene_whole_W200N50F50_chr", 
    data_type = "rpkm/mean"))
```

```
## [1] "/media/storage2/analysis/profiles/norm/rpkm/mean/gene_whole_W200N50F50_chr/images/omp_hmc_rep1_q30_rmdup_extend300_mean_omp_hmc_rep2_q30_rmdup"
## [1] "/media/storage2/analysis/profiles/norm/rpkm/mean/gene_whole_W200N50F50_chr/images/d3xog_het_hmc_sort_q30_rmdup"
## [1] "/media/storage2/analysis/profiles/norm/rpkm/mean/gene_whole_W200N50F50_chr/images/d3xog_ko_hmc_sort_q30_rmdup"
## [1] "/media/storage2/analysis/profiles/norm/rpkm/mean/gene_whole_W200N50F50_chr/images/omp_mc_rep1_q30_rmdup_extend300"
## [1] "/media/storage2/analysis/profiles/norm/rpkm/mean/gene_whole_W200N50F50_chr/images/d3xog_het_mc_sort_q30_rmdup"
## [1] "/media/storage2/analysis/profiles/norm/rpkm/mean/gene_whole_W200N50F50_chr/images/d3xog_ko_mc_sort_q30_rmdup"
## [1] "/media/storage2/analysis/profiles/norm/rpkm/mean/gene_whole_W200N50F50_chr/images/ngn_hmc_rep1_q30_rmdup_extend300_mean_ngn_hmc_rep2_q30_rmdup"
## [1] "/media/storage2/analysis/profiles/norm/rpkm/mean/gene_whole_W200N50F50_chr/images/icam_hmc_rep1_q30_rmdup_extend300_mean_icam_hmc_rep2_q30_rmdup"
## [1] "/media/storage2/analysis/profiles/norm/rpkm/mean/gene_whole_W200N50F50_chr/images/ngn_mc_rep1_q30_rmdup_extend300_mean_ngn_mc_rep2_q30_rmdup"
## [1] "/media/storage2/analysis/profiles/norm/rpkm/mean/gene_whole_W200N50F50_chr/images/icam_mc_rep1_q30_rmdup_extend300_mean_icam_mc_rep2_q30_rmdup"
```


#### Compute means of gene body levels

```r
data.mid <- lapply(data, function(x) apply(x[, 51:100], 1, mean))
data.mid <- data.frame(do.call("cbind", data.mid))
colnames(data.mid) <- c("wt_omp_hmc", "het_omp_hmc", "ko_omp_hmc", "wt_omp_mc", 
    "het_omp_mc", "ko_omp_mc", "wt_ngn_hmc", "wt_icam_hmc", "wt_ngn_mc", "wt_icam_mc")
data.mid <- data.mid[order(data.mid[, 1]), ]
```



```r
data.mid$id <- rownames(data.mid)
data.mid.m <- melt(data.mid)
s <- str_split(data.mid.m$variable, "_")
data.mid.m$geno <- factor(unlist(lapply(s, function(x) x[1])), levels = c("wt", 
    "het", "ko"))
levels(data.mid.m$geno) <- c("+/+", "+/-", "-/-")
data.mid.m$celltype <- factor(unlist(lapply(s, function(x) x[2])), levels = c("omp", 
    "ngn", "icam"))
levels(data.mid.m$celltype) <- c("mOSN", "GBC", "HBC")
data.mid.m$mod <- factor(unlist(lapply(s, function(x) x[3])), levels = c("hmc", 
    "mc"))
levels(data.mid.m$mod) <- c("5hmC", "5mC")
```



```r
gene <- read.delim("~/s2/data/rna/cuffdiff/d3xog_wt_het_ko_rmrna_run9/gene_exp.diff")
gene.read <- read.delim("~/s2/data/rna/cuffdiff/d3xog_wt_het_ko_rmrna_run9/genes.read_group_tracking")
gene.read.fpkm <- dcast(gene.read, tracking_id ~ condition + replicate, value.var = "FPKM")
gene.read.fpkm$gene <- gene$gene[match(gene.read.fpkm$tracking_id, gene$test_id)]
gene.read.fpkm <- na.omit(gene.read.fpkm[match(rownames(data.mid), gene.read.fpkm$gene), 
    ])
gene.read.fpkm <- numcolwise(onelog2)(gene.read.fpkm)
cor(gene.read.fpkm)
```

```
##        q1_0   q1_1   q1_2   q2_0   q3_0   q3_1
## q1_0 1.0000 0.9721 0.9436 0.9752 0.9178 0.9589
## q1_1 0.9721 1.0000 0.9531 0.9669 0.9097 0.9721
## q1_2 0.9436 0.9531 1.0000 0.9462 0.9377 0.9555
## q2_0 0.9752 0.9669 0.9462 1.0000 0.9264 0.9645
## q3_0 0.9178 0.9097 0.9377 0.9264 1.0000 0.9186
## q3_1 0.9589 0.9721 0.9555 0.9645 0.9186 1.0000
```

```r
gene.read.fpkm.m <- melt(gene.read.fpkm)
levels(gene.read.fpkm.m$variable) <- c("+/+", "+/-", "-/- rep1", "-/- rep2")
```



```r
gg <- ggplot(gene.read.fpkm.m, aes(value, color = variable))
gg + geom_density() + facet_grid(variable ~ .) + scale_color_brewer(palette = "Set1") + 
    labs(x = "log2(FPKM+1)")
```

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6.png) 



```r
gene.read.fpkm.0 <- gene.read.fpkm[, c("q3_0", "q2_0", "q1_0")]
gene.read.fpkm.0.sub <- pairwise(gene.read.fpkm.0, subtract)
```

```
## [1] "q3_0" "q2_0"
## [1] "q3_0" "q1_0"
## [1] "q2_0" "q1_0"
```



```r
gene.read.fpkm.c100 <- chunkMatrix(gene.read.fpkm, 100, median)
head(gene.read.fpkm.c100)
```

```
##          q1_0 q1_1 q1_2 q2_0   q3_0   q3_1 index
## result.1    0    0    0    0 0.0000 0.0000     1
## result.2    0    0    0    0 0.0000 0.0000     2
## result.3    0    0    0    0 0.0000 0.0000     3
## result.4    0    0    0    0 0.0000 0.1481     4
## result.5    0    0    0    0 0.0000 0.1047     5
## result.6    0    0    0    0 0.1329 0.1373     6
```

```r
gene.read.fpkm.c100.m <- melt(gene.read.fpkm.c100, id.vars = "index")
levels(gene.read.fpkm.c100.m$variable) <- c("+/+", "+/-", "-/- rep1", "-/- rep2")
```


### Replicates separate
#### Ordered by gene body 5hmC 




































































