D3a - DNase peaks 2 - heatmaps
========================================================


```r
# opts_chunk$set(warning=FALSE, message=FALSE, error=FALSE)
suppressPackageStartupMessages(source("~/src/seqAnalysis/R/profiles2.R"))
suppressPackageStartupMessages(source("~/src/seqAnalysis/R/image.R"))
library(ggplot2)
library(reshape2)
library(gridExtra)
```



```r
positionMatrix.all("d3a_het_dnase_sort_q30_dnase_sub_igenome_ensembl_genes_extend5kb_nosex.bed_W25F200_both_chr", 
    data_type = "rpkm/mean")
```



```r
wt.hmc <- makeImage("omp_hmc_rep1_q30_rmdup_extend300_mean_omp_hmc_rep2_q30_rmdup", 
    "d3a_het_dnase_sort_q30_dnase_sub_igenome_ensembl_genes_extend5kb_nosex.bed_W25F200_both_chr", 
    data_type = "rpkm/mean", image = F)
```

```
## [1] "/media/storage2/analysis/profiles/norm/rpkm/mean/d3a_het_dnase_sort_q30_dnase_sub_igenome_ensembl_genes_extend5kb_nosex.bed_W25F200_both_chr/images/omp_hmc_rep1_q30_rmdup_extend300_mean_omp_hmc_rep2_q30_rmdup"
```

```r
ko.hmc <- makeImage("d3xog_ko_hmc_paired_q30", "d3a_het_dnase_sort_q30_dnase_sub_igenome_ensembl_genes_extend5kb_nosex.bed_W25F200_both_chr", 
    data_type = "rpkm/mean", image = F)
```

```
## [1] "/media/storage2/analysis/profiles/norm/rpkm/mean/d3a_het_dnase_sort_q30_dnase_sub_igenome_ensembl_genes_extend5kb_nosex.bed_W25F200_both_chr/images/d3xog_ko_hmc_paired_q30"
```

```r

n <- which(is.na(wt.hmc), arr.ind = T)
wt.hmc[n] <- 0
n <- which(is.na(ko.hmc), arr.ind = T)
ko.hmc[n] <- 0

ko.wt.hmc <- ko.hmc - wt.hmc
```



```r
wt.hmc.pc <- prcomp(wt.hmc)
wt.hmc.pred <- predict(wt.hmc.pc, wt.hmc)
```



```r
MP.heat(wt.hmc[order(wt.hmc.pred[, 1]), ], average = 20, range = c(0, 3))
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5.png) 



```r
MP.heat(ko.hmc[order(wt.hmc.pred[, 1]), ], average = 20, range = c(0, 3))
```

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6.png) 


Middle 500 bp 

```r
wt.hmc.mid <- apply(wt.hmc[, 191:210], 1, mean)
ko.hmc.mid <- apply(ko.hmc[, 191:210], 1, mean)

mid <- cbind(wt = wt.hmc.mid, ko = ko.hmc.mid)
mid <- mid[order(mid[, 1]), ]
mid.c100 <- chunkMatrix(mid)
```

```
## Error: could not find function "chunkMatrix"
```

```r
mid.c100.m <- melt(mid.c100, id.vars = "index")
```

```
## Error: object 'mid.c100' not found
```



```r
gg <- ggplot(mid.c100.m, aes(index, value))
```

```
## Error: object 'mid.c100.m' not found
```

```r
gg.mid <- gg + geom_point(aes(color = variable)) + scale_color_manual(values = col2) + 
    labs(title = "Middle")
```

```
## Error: object 'gg' not found
```


Flanking (-1kb to -500, +500 to +1kb)

```r
wt.hmc.flank <- apply(wt.hmc[, c(171:190, 211:230)], 1, mean)
ko.hmc.flank <- apply(ko.hmc[, c(171:190, 211:230)], 1, mean)

flank <- cbind(wt = wt.hmc.flank, ko = ko.hmc.flank)
flank <- flank[order(flank[, 1]), ]
flank.c100 <- chunkMatrix(flank)
```

```
## Error: could not find function "chunkMatrix"
```

```r
flank.c100.m <- melt(flank.c100, id.vars = "index")
```

```
## Error: object 'flank.c100' not found
```



```r
gg <- ggplot(flank.c100.m, aes(index, value))
```

```
## Error: object 'flank.c100.m' not found
```

```r
gg.flank <- gg + geom_point(aes(color = variable)) + scale_color_manual(values = col2) + 
    labs(title = "Flanking")
```

```
## Error: object 'gg' not found
```



```r
grid.arrange(gg.mid, gg.flank, ncol = 2)
```

```
## Error: object 'gg.mid' not found
```



```r
wt.nuc <- makeImage("d3xog_wt_nuc_478_rmdup", "d3a_het_dnase_sort_q30_dnase_sub_igenome_ensembl_genes_extend5kb_nosex.bed_W25F200_both_chr", 
    data_type = "rpkm/mean", image = F)
```

```
## [1] "/media/storage2/analysis/profiles/norm/rpkm/mean/d3a_het_dnase_sort_q30_dnase_sub_igenome_ensembl_genes_extend5kb_nosex.bed_W25F200_both_chr/images/d3xog_wt_nuc_478_rmdup"
```

```r
ko.nuc <- makeImage("d3xog_ko_nuc_256_rmdup", "d3a_het_dnase_sort_q30_dnase_sub_igenome_ensembl_genes_extend5kb_nosex.bed_W25F200_both_chr", 
    data_type = "rpkm/mean", image = F)
```

```
## [1] "/media/storage2/analysis/profiles/norm/rpkm/mean/d3a_het_dnase_sort_q30_dnase_sub_igenome_ensembl_genes_extend5kb_nosex.bed_W25F200_both_chr/images/d3xog_ko_nuc_256_rmdup"
```

```r

n <- which(is.na(wt.nuc), arr.ind = T)
wt.nuc[n] <- 0
n <- which(is.na(ko.nuc), arr.ind = T)
ko.nuc[n] <- 0

ko.wt.nuc <- ko.nuc - wt.nuc
```



```r
MP.heat(wt.nuc[order(wt.hmc.pred[, 1]), ], average = 20, range = c(0, 1))
```

![plot of chunk unnamed-chunk-13](figure/unnamed-chunk-13.png) 



```r
MP.heat(ko.nuc[order(wt.hmc.pred[, 1]), ], average = 20, range = c(0, 1))
```

![plot of chunk unnamed-chunk-14](figure/unnamed-chunk-14.png) 

Middle 100 bp 

```r
wt.nuc.mid <- apply(wt.nuc[, 193:204], 1, mean)
ko.nuc.mid <- apply(ko.nuc[, 193:204], 1, mean)

mid <- cbind(wt = wt.nuc.mid, ko = ko.nuc.mid)
mid <- mid[order(mid[, 1]), ]
mid.c100 <- chunkMatrix(mid)
```

```
## Error: could not find function "chunkMatrix"
```

```r
mid.c100.m <- melt(mid.c100, id.vars = "index")
```

```
## Error: object 'mid.c100' not found
```



```r
gg <- ggplot(mid.c100.m, aes(index, value))
```

```
## Error: object 'mid.c100.m' not found
```

```r
gg.mid1 <- gg + geom_point(aes(color = variable)) + scale_color_manual(values = col2) + 
    labs(title = "Middle") + coord_cartesian(ylim = c(0, 1))
```

```
## Error: object 'gg' not found
```



Middle 500 bp 

```r
wt.nuc.mid <- apply(wt.nuc[, 191:210], 1, mean)
ko.nuc.mid <- apply(ko.nuc[, 191:210], 1, mean)

mid <- cbind(wt = wt.nuc.mid, ko = ko.nuc.mid)
mid <- mid[order(mid[, 1]), ]
mid.c100 <- chunkMatrix(mid)
```

```
## Error: could not find function "chunkMatrix"
```

```r
mid.c100.m <- melt(mid.c100, id.vars = "index")
```

```
## Error: object 'mid.c100' not found
```



```r
gg <- ggplot(mid.c100.m, aes(index, value))
```

```
## Error: object 'mid.c100.m' not found
```

```r
gg.mid <- gg + geom_point(aes(color = variable)) + scale_color_manual(values = col2) + 
    labs(title = "Middle") + coord_cartesian(ylim = c(0, 1))
```

```
## Error: object 'gg' not found
```


Flanking (-1kb to -500, +500 to +1kb)

```r
wt.nuc.flank <- apply(wt.nuc[, c(171:190, 211:230)], 1, mean)
ko.nuc.flank <- apply(ko.nuc[, c(171:190, 211:230)], 1, mean)

flank <- cbind(wt = wt.nuc.flank, ko = ko.nuc.flank)
flank <- flank[order(flank[, 1]), ]
flank.c100 <- chunkMatrix(flank)
```

```
## Error: could not find function "chunkMatrix"
```

```r
flank.c100.m <- melt(flank.c100, id.vars = "index")
```

```
## Error: object 'flank.c100' not found
```



```r
gg <- ggplot(flank.c100.m, aes(index, value))
```

```
## Error: object 'flank.c100.m' not found
```

```r
gg.flank <- gg + geom_point(aes(color = variable)) + scale_color_manual(values = col2) + 
    labs(title = "Flanking") + coord_cartesian(ylim = c(0, 1))
```

```
## Error: object 'gg' not found
```



```r
grid.arrange(gg.mid1, gg.mid, gg.flank, ncol = 3)
```

```
## Error: object 'gg.mid1' not found
```

