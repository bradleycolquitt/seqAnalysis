D3a - DNase peaks 2 - heatmaps
========================================================


```r
suppressPackageStartupMessages(source("~/src/seqAnalysis/R/profiles2.R"))
suppressPackageStartupMessages(source("~/src/seqAnalysis/R/image.R"))
suppressPackageStartupMessages(source("~/src/seqAnalysis/R/features.R"))
library(ggplot2)
library(reshape2)
```

```
## Attaching package: 'reshape2'
```

```
## The following object(s) are masked from 'package:reshape':
## 
## colsplit, melt, recast
```

```r
library(gridExtra)
```


### Peak

```r
cells <- makeFeatureMatrix2("d3a_het_dnase_inter_moe_h3k27ac_interV_moe_h3k4me.bed_chr", 
    "cells_rep_hmc", data_type = "rpkm/mean")
```

```
## [1] "omp_hmc_rep1_q30_rmdup_extend300_mean_omp_hmc_rep2_q30_rmdup"  
## [2] "ngn_hmc_rep1_q30_rmdup_extend300_mean_ngn_hmc_rep2_q30_rmdup"  
## [3] "icam_hmc_rep1_q30_rmdup_extend300_mean_icam_hmc_rep2_q30_rmdup"
```

```r
d3a <- makeFeatureMatrix2("d3a_het_dnase_inter_moe_h3k27ac_interV_moe_h3k4me.bed_chr", 
    "d3xog_hmc", data_type = "rpkm/mean")
```

```
## [1] "omp_hmc_rep1_q30_rmdup_extend300_mean_omp_hmc_rep2_q30_rmdup"
## [2] "d3xog_het_hmc_paired_q30"                                    
## [3] "d3xog_ko_hmc_paired_q30"
```



```r
comb <- cbind(cells, d3a[, 3])
```



```r
comb <- comb[order(comb[, "omp_hmc"]), ]
comb.c100 <- chunkMatrix(comb, chunks = 100)
comb.c100.m <- melt(comb.c100, id.vars = "index")
levels(comb.c100.m$variable) <- c("mOSN WT", "GBC", "HBC", "mOSN KO")
```



```r
theme_set(theme_classic())
gg <- ggplot(comb.c100.m, aes(index, value, color = variable))
gg + geom_point() + scale_color_manual(values = c("skyblue4", "slategray4", 
    "slategray3", "tomato4")) + labs(y = "Average RPM") + theme(legend.position = c(0.2, 
    0.8), legend.title = element_blank(), legend.key = element_blank())
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5.png) 


### Flank +- 1kb

```r
cells <- makeFeatureMatrix2("d3a_het_dnase_inter_moe_h3k27ac_interV_moe_h3k4me.bed_chr", 
    "cells_rep_hmc", data_type = "rpkm/mean", select = "flank1000")
```

```
##               omp_hmc_rep1_q30_rmdup_extend300_mean_omp_hmc_rep2_q30_rmdup 
##   "omp_hmc_rep1_q30_rmdup_extend300_mean_omp_hmc_rep2_q30_rmdup_flank1000" 
##               ngn_hmc_rep1_q30_rmdup_extend300_mean_ngn_hmc_rep2_q30_rmdup 
##   "ngn_hmc_rep1_q30_rmdup_extend300_mean_ngn_hmc_rep2_q30_rmdup_flank1000" 
##             icam_hmc_rep1_q30_rmdup_extend300_mean_icam_hmc_rep2_q30_rmdup 
## "icam_hmc_rep1_q30_rmdup_extend300_mean_icam_hmc_rep2_q30_rmdup_flank1000"
```

```r
d3a <- makeFeatureMatrix2("d3a_het_dnase_inter_moe_h3k27ac_interV_moe_h3k4me.bed_chr", 
    "d3xog_hmc", data_type = "rpkm/mean", select = "flank1000")
```

```
##             omp_hmc_rep1_q30_rmdup_extend300_mean_omp_hmc_rep2_q30_rmdup 
## "omp_hmc_rep1_q30_rmdup_extend300_mean_omp_hmc_rep2_q30_rmdup_flank1000" 
##                                                 d3xog_het_hmc_paired_q30 
##                                     "d3xog_het_hmc_paired_q30_flank1000" 
##                                                  d3xog_ko_hmc_paired_q30 
##                                      "d3xog_ko_hmc_paired_q30_flank1000"
```



```r
comb <- cbind(cells, d3a[, 3])
```



```r
comb <- comb[order(comb[, "omp_hmc"]), ]
comb.c100 <- chunkMatrix(comb, chunks = 100)
comb.c100.m <- melt(comb.c100, id.vars = "index")
levels(comb.c100.m$variable) <- c("mOSN WT", "GBC", "HBC", "mOSN KO")
```



```r
theme_set(theme_classic())
gg <- ggplot(comb.c100.m, aes(index, value, color = variable))
gg + geom_point() + scale_color_manual(values = c("skyblue4", "slategray4", 
    "slategray3", "tomato4")) + labs(y = "Average RPM") + theme(legend.position = c(0.2, 
    0.8), legend.title = element_blank(), legend.key = element_blank())
```

![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-9.png) 



```r
grid.arrange(gg.mid1, gg.mid, gg.flank, ncol = 3)
```

```
## Error: object 'gg.mid1' not found
```

