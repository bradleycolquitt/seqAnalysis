D3xOG rmrna -- run6
========================================================

Cuffdiff parameters:
  * version 2.1.1
  * -b /seq/lib/indexes/mm9.fa 
  * -M /home/user/lib/rmsk/rmsk.gtf 
  * --compatible-hits-norm
  * --multi-read-correct
  * library-norm-method classic-fpkm
  * dispersion-method pooled
  * genes_chr_protein.gtf

Samples
  * wt
    * 121126/omp_rmrna_blank/omp_rmrna_rep1_blank.bam
    * 130326/omp_rmrna_rep2_blank/omp_rmrna_rep2_blank.bam
    * 121126/d3xog_wt_rmrna_blank/d3xog_wt_rmrna_blank.bam
  * ko
    * 121126/d3xog_ko_rmrna_blank/d3xog_ko_rmrna_blank.bam
    * 130326/d3xog_ko_rmrna_rep2_blank/d3xog_ko_rmrna_rep2_blank.bam
    

```r
opts_chunk$set(warning = FALSE, message = FALSE, error = FALSE)
library(plyr)
library(reshape2)
library(ggplot2)
suppressPackageStartupMessages(source("~/src/seqAnalysis/R/features.R"))
suppressPackageStartupMessages(source("~/src/seqAnalysis/R/image.R"))
```



```r
wt.hmc <- makeImage("omp_hmc_rep1_mean_omp_hmc_rep2", "refGene_noRandom_order_outsides2_tss_W25F200_chr", 
    data_type = "rpkm/mean", image = F)
```

```
## [1] "/media/storage2/analysis/profiles/norm/rpkm/mean/refGene_noRandom_order_outsides2_tss_W25F200_chr/images/omp_hmc_rep1_mean_omp_hmc_rep2"
```

```r
# wt.hmc <- makeImage('omp_hmc_120424_rpkm',
# 'refGene_noRandom_order_outsides2_tss_W25F200_chr',
# data_type='rpkm/mean', image=F)
het.hmc <- makeImage("d3xog_het_hmc", "refGene_noRandom_order_outsides2_tss_W25F200_chr", 
    data_type = "rpkm/mean", image = F)
```

```
## [1] "/media/storage2/analysis/profiles/norm/rpkm/mean/refGene_noRandom_order_outsides2_tss_W25F200_chr/images/d3xog_het_hmc"
```

```r
ko.hmc <- makeImage("d3xog_ko_hmc", "refGene_noRandom_order_outsides2_tss_W25F200_chr", 
    data_type = "rpkm/mean", image = F)
```

```
## [1] "/media/storage2/analysis/profiles/norm/rpkm/mean/refGene_noRandom_order_outsides2_tss_W25F200_chr/images/d3xog_ko_hmc"
```


#### Compute means of region +1kb to +2kb of TSS

```r
wt.hmc.up <- apply(wt.hmc[, 241:280], 1, mean)
het.hmc.up <- apply(het.hmc[, 241:280], 1, mean)
ko.hmc.up <- apply(ko.hmc[, 241:280], 1, mean)

wt.hmc.up <- wt.hmc.up[order(wt.hmc.up)]
het.hmc.up <- het.hmc.up[match(names(wt.hmc.up), names(het.hmc.up))]
ko.hmc.up <- ko.hmc.up[match(names(wt.hmc.up), names(ko.hmc.up))]

hmc.up <- cbind(wt.hmc.up, het.hmc.up, ko.hmc.up)
hmc.up.c100 <- chunkMatrix(hmc.up, 100)
hmc.up.c100.m <- melt(hmc.up.c100, id.vars = "index")
```



```r
gene <- read.delim("~/s2/data/rna/cuffdiff/d3xog_wt_het_ko_rmrna_run6/gene_exp.diff")
gene.read <- read.delim("~/s2/data/rna/cuffdiff/d3xog_wt_het_ko_rmrna_run6/genes.read_group_tracking")
gene.read.fpkm <- dcast(gene.read, tracking_id ~ condition + replicate, value.var = "FPKM")
gene.read.fpkm$gene <- gene$gene[match(gene.read.fpkm$tracking_id, gene$test_id)]
gene.read.fpkm <- na.omit(gene.read.fpkm[match(names(wt.hmc.up), gene.read.fpkm$gene), 
    ])
gene.read.fpkm <- numcolwise(onelog2)(gene.read.fpkm)
```



```r
par(mfrow = c(2, 3))
apply(gene.read.fpkm, 2, function(x) plot(density(x)))
```

```
## NULL
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5.png) 



```r
gene.read.fpkm.c100 <- chunkMatrix(gene.read.fpkm, 100)
gene.read.fpkm.c100.m <- melt(gene.read.fpkm.c100, id.vars = "index")
```


#### Replicates separate

```r
gg.rep <- ggplot(gene.read.fpkm.c100.m, aes(index, value, color = variable))
gg.rep + geom_point() + scale_color_manual(values = c("red1", "red2", "red3", 
    "green1", "green2")) + facet_wrap(~variable)
```

![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-7.png) 



```r
gene.fpkm <- read.delim("~/s2/data/rna/cuffdiff/d3xog_wt_het_ko_rmrna_run6/genes.fpkm_tracking")
gene.fpkm.mat <- gene.fpkm[, c("gene_short_name", "q1_FPKM", "q2_FPKM")]
gene.fpkm.mat <- na.omit(gene.fpkm.mat[match(names(wt.hmc.up), gene.fpkm.mat[, 
    1]), ])
gene.fpkm.mat <- numcolwise(onelog2)(gene.fpkm.mat)
gene.fpkm.mat <- transform(gene.fpkm.mat, q2.q1 = q2_FPKM - q1_FPKM)
gene.fpkm.mat.c100 <- chunkMatrix(gene.fpkm.mat, 100)
gene.fpkm.mat.c100.m <- melt(gene.fpkm.mat.c100, id.vars = "index")
```



```r
gg <- ggplot(gene.fpkm.mat.c100.m, aes(index, value, color = variable))
gg + geom_point()
```

![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-9.png) 



```r
gene.fpkm.mat$q2.q1 <- with(gene.fpkm.mat, q2_FPKM - q1_FPKM)
par(mfrow = c(1, 3))
apply(gene.fpkm.mat, 2, function(x) plot(density(x)))
```

![plot of chunk unnamed-chunk-10](figure/unnamed-chunk-10.png) 

```
## NULL
```

