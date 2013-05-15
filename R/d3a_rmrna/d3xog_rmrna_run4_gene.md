D3xOG rmrna -- run4
========================================================

Cuffdiff parameters:
  * version 2.1.1
  * -b /seq/lib/indexes/mm9.fa 
  * -M /home/user/lib/rmsk/rmsk.gtf 
  * --compatible-hits-norm
  * -u
  * --library-norm-method geometric
  * dispersion-method per-condition
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
wt.hmc <- makeImage("omp_hmc_rep1_mean_omp_hmc_rep2", "gene_whole_W200N50F50_chr", 
    data_type = "rpkm/mean", image = F)
```

```
## [1] "/media/storage2/analysis/profiles/norm/rpkm/mean/gene_whole_W200N50F50_chr/images/omp_hmc_rep1_mean_omp_hmc_rep2"
```

```r
# wt.hmc <- makeImage('omp_hmc_120424_rpkm', 'gene_whole_W200N50F50_chr',
# data_type='rpkm/mean', image=F)
het.hmc <- makeImage("d3xog_het_hmc_paired_q30", "gene_whole_W200N50F50_chr", 
    data_type = "rpkm/mean", image = F)
```

```
## [1] "/media/storage2/analysis/profiles/norm/rpkm/mean/gene_whole_W200N50F50_chr/images/d3xog_het_hmc_paired_q30"
```

```r
ko.hmc <- makeImage("d3xog_ko_hmc_paired_q30", "gene_whole_W200N50F50_chr", 
    data_type = "rpkm/mean", image = F)
```

```
## [1] "/media/storage2/analysis/profiles/norm/rpkm/mean/gene_whole_W200N50F50_chr/images/d3xog_ko_hmc_paired_q30"
```


#### Compute means of gene body levels

```r
wt.hmc.up <- apply(wt.hmc[, 51:100], 1, mean)
het.hmc.up <- apply(het.hmc[, 51:100], 1, mean)
ko.hmc.up <- apply(ko.hmc[, 51:100], 1, mean)


wt.hmc.up <- wt.hmc.up[order(wt.hmc.up)]
het.hmc.up <- het.hmc.up[match(names(wt.hmc.up), names(het.hmc.up))]
ko.hmc.up <- ko.hmc.up[match(names(wt.hmc.up), names(ko.hmc.up))]

hmc.up <- cbind(wt.hmc.up, het.hmc.up, ko.hmc.up)
hmc.up.c100 <- chunkMatrix(hmc.up, 100)
hmc.up.c100.m <- melt(hmc.up.c100, id.vars = "index")
```



```r
gg <- ggplot(hmc.up.c100.m, aes(index, value, color = variable))
gg + geom_point() + coord_cartesian(xlim = c(20, 100))
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4.png) 



```r
gene <- read.delim("~/s2/data/rna/cuffdiff/d3xog_wt_het_ko_rmrna_run4/gene_exp.diff")
gene.read <- read.delim("~/s2/data/rna/cuffdiff/d3xog_wt_het_ko_rmrna_run4/genes.read_group_tracking")
gene.read.fpkm <- dcast(gene.read, tracking_id ~ condition + replicate, value.var = "FPKM")
gene.read.fpkm$gene <- gene$gene[match(gene.read.fpkm$tracking_id, gene$test_id)]
gene.read.fpkm <- na.omit(gene.read.fpkm[match(names(wt.hmc.up), gene.read.fpkm$gene), 
    ])
gene.read.fpkm <- numcolwise(onelog2)(gene.read.fpkm)
```



```r
gene.read.fpkm.c100 <- chunkMatrix(gene.read.fpkm, 100)
gene.read.fpkm.c100.m <- melt(gene.read.fpkm.c100, id.vars = "index")
```



```r
gene.read.fpkm.0 <- gene.read.fpkm[, c("q1_0", "q2_0")]
gene.read.fpkm.0.sub <- subtract(gene.read.fpkm.0[, 2], gene.read.fpkm.0[, 1])
gene.read.fpkm.0$q2_0_q1_0 <- gene.read.fpkm.0.sub
```



```r
gene.read.fpkm.1 <- gene.read.fpkm[, c("q1_0", "q2_1")]
gene.read.fpkm.1.sub <- subtract(gene.read.fpkm.1[, 2], gene.read.fpkm.1[, 1])
gene.read.fpkm.1$q2_1_q1_0 <- gene.read.fpkm.1.sub
```


#### Replicates separate

```r
gg.rep <- ggplot(gene.read.fpkm.c100.m, aes(index, value, color = variable))
gg.rep + geom_point() + scale_color_manual(values = c("red1", "red2", "red3", 
    "green1", "green2")) + facet_wrap(~variable)
```

![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-9.png) 



```r
gene.fpkm <- read.delim("~/s2/data/rna/cuffdiff/d3xog_wt_het_ko_rmrna_run4/genes.fpkm_tracking")
gene.fpkm.mat <- gene.fpkm[, c("gene_short_name", "q1_FPKM", "q2_FPKM")]
gene.fpkm.mat <- na.omit(gene.fpkm.mat[match(names(wt.hmc.up), gene.fpkm.mat[, 
    1]), ])
gene.fpkm.mat <- numcolwise(onelog2)(gene.fpkm.mat)
gene.fpkm.mat <- transform(gene.fpkm.mat, q2.q1 = q2_FPKM - q1_FPKM)
gene.fpkm.mat.c100 <- chunkMatrix(gene.fpkm.mat, 100, median)
gene.fpkm.mat.c100.m <- melt(gene.fpkm.mat.c100, id.vars = "index")
```




```r
gg <- ggplot(gene.fpkm.mat.c100.m, aes(index, value, color = variable))
gg + geom_point()
```

![plot of chunk unnamed-chunk-11](figure/unnamed-chunk-11.png) 



```r
gene.read.fpkm.0.c100 <- chunkMatrix(gene.read.fpkm.0, 100, mean)
head(gene.read.fpkm.0.c100)
```

```
##             q1_0    q2_0 q2_0_q1_0 index
## result.1 0.08195 0.07997 -0.001976     1
## result.2 0.16101 0.20889  0.047887     2
## result.3 0.16159 0.24710  0.085508     3
## result.4 0.25120 0.50090  0.249704     4
## result.5 0.23256 0.43453  0.201962     5
## result.6 0.47203 0.59842  0.126381     6
```

```r
gene.read.fpkm.0.c100.m <- melt(gene.read.fpkm.0.c100, id.vars = "index")
```



```r
gg.rep <- ggplot(gene.read.fpkm.0.c100.m, aes(index, value, color = variable))
gg.rep + geom_point() + scale_color_manual(values = c("red1", "green2", "purple"))
```

![plot of chunk unnamed-chunk-13](figure/unnamed-chunk-13.png) 



```r
gene.read.fpkm.1.c100 <- chunkMatrix(gene.read.fpkm.1, 100, mean)
head(gene.read.fpkm.1.c100)
```

```
##             q1_0    q2_1 q2_1_q1_0 index
## result.1 0.08195 0.08891  0.006959     1
## result.2 0.16101 0.16226  0.001254     2
## result.3 0.16159 0.14820 -0.013393     3
## result.4 0.25120 0.36916  0.117960     4
## result.5 0.23256 0.37886  0.146298     5
## result.6 0.47203 0.53869  0.066659     6
```

```r
gene.read.fpkm.1.c100.m <- melt(gene.read.fpkm.1.c100, id.vars = "index")
```



```r
gg.rep <- ggplot(gene.read.fpkm.1.c100.m, aes(index, value, color = variable))
gg.rep + geom_point() + scale_color_manual(values = c("red1", "green2", "purple"))
```

![plot of chunk unnamed-chunk-15](figure/unnamed-chunk-15.png) 



```r
par(mfrow = c(1, 3))
apply(gene.fpkm.mat, 2, function(x) plot(density(x)))
```

![plot of chunk unnamed-chunk-16](figure/unnamed-chunk-16.png) 

```
## NULL
```



