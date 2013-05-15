D3xOG 5hmC/5mC TSS
========================================================


```r
opts_chunk$set(warning = FALSE, message = FALSE, error = FALSE)
library(plyr)
library(reshape2)
library(ggplot2)
suppressPackageStartupMessages(source("~/src/seqAnalysis/R/profiles2.R"))
suppressPackageStartupMessages(source("~/src/seqAnalysis/R/image.R"))
```



```r
makeProfile2.allSamp("refGene_noRandom_order_outsides2_tss_W25F200_chr", data_type = "rpkm/mean", 
    rm.outliers = 0.01)
```

```
## [1] "/media/storage2/analysis/profiles/norm/rpkm/mean/refGene_noRandom_order_outsides2_tss_W25F200_chr"
## Note: next may be used in wrong context: no loop is visible
```



```r
plot2.several("refGene_noRandom_order_outsides2_tss_W25F200_chr", "d3xog_hmc", 
    data_type = "rpkm/mean", cols = col3, group2 = "trim0.01")
```

```
## [1] "omp_hmc_rep1_mean_omp_hmc_rep2_trim0.01"
## [1] "omp_hmc_rep1_mean_omp_hmc_rep2_trim0.01_mean"
## [1] "d3xog_het_hmc_paired_q30_trim0.01"
## [1] "d3xog_het_hmc_paired_q30_trim0.01_mean"
## [1] "d3xog_ko_hmc_paired_q30_trim0.01"
## [1] "d3xog_ko_hmc_paired_q30_trim0.01_mean"
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3.png) 

```
## [1] 0.093 0.820
```



```r
plot2.several("refGene_noRandom_order_outsides2_tss_W25F200_chr", "d3xog_mc", 
    data_type = "rpkm/mean", cols = col3, group2 = "trim0.01")
```

```
## [1] "omp_mc_rmdup_trim0.01"
## [1] "omp_mc_rmdup_trim0.01_mean"
## [1] "d3xog_het_mc_paired_q30_trim0.01"
## [1] "d3xog_het_mc_paired_q30_trim0.01_mean"
## [1] "d3xog_ko_mc_paired_q30_trim0.01"
## [1] "d3xog_ko_mc_paired_q30_trim0.01_mean"
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4.png) 

```
## [1] 0.136 0.574
```



```r
plot2.several("refGene_noRandom_order_outsides2_tss_W25F200_chr", "tt3_rep", 
    data_type = "rpkm/mean", cols = col2, group2 = "trim0.01")
```

```
## [1] "omp_hmc_rep1_mean_omp_hmc_rep2_trim0.01"
## [1] "omp_hmc_rep1_mean_omp_hmc_rep2_trim0.01_mean"
## [1] "ott3_hmc_rep1_mean_ott3_hmc_rep2_trim0.01"
## [1] "ott3_hmc_rep1_mean_ott3_hmc_rep2_trim0.01_mean"
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5.png) 

```
## [1] 0.048 0.829
```



```r
makeProfile2.allSamp("gene_whole_W200N50F50_chr", data_type = "rpkm/mean", rm.outliers = 0.01)
```

```
## [1] "/media/storage2/analysis/profiles/norm/rpkm/mean/gene_whole_W200N50F50_chr"
## Note: next may be used in wrong context: no loop is visible
```



```r
plot2.several("gene_whole_W200N50F50_chr", "d3a_hmc", data_type = "rpkm/mean", 
    cols = col3[c(1, 3)], group2 = "trim0.01", y.vals = c(0, 0.8))
```

```
## [1] "moe_d3a_wt_hmc_rpkm_trim0.01"
## [1] "moe_d3a_wt_hmc_rpkm_trim0.01_mean"
## [1] "moe_d3a_ko_hmc_rpkm_trim0.01"
## [1] "moe_d3a_ko_hmc_rpkm_trim0.01_mean"
```

![plot of chunk gene_whole_W200N50F50_d3a_hmc](figure/gene_whole_W200N50F50_d3a_hmc.png) 

```
## [1] 0.0 0.8
```



```r
plot2.several("gene_whole_W200N50F50_chr", "d3a_mc", data_type = "rpkm/mean", 
    cols = col3[c(1, 3)], group2 = "trim0.01", y.vals = c(0, 0.8))
```

```
## [1] "moe_d3a_wt_mc_rpkm_trim0.01"
## [1] "moe_d3a_wt_mc_rpkm_trim0.01_mean"
## [1] "moe_d3a_ko_mc_rpkm_trim0.01"
## [1] "moe_d3a_ko_mc_rpkm_trim0.01_mean"
```

![plot of chunk gene_whole_W200N50F50_d3a_mc](figure/gene_whole_W200N50F50_d3a_mc.png) 

```
## [1] 0.0 0.8
```



```r
plot2.several("gene_whole_W200N50F50_chr", "d3xog_hmc", data_type = "rpkm/mean", 
    cols = col3, group2 = "trim0.01", y.vals = c(0, 0.8))
```

```
## [1] "omp_hmc_rep1_mean_omp_hmc_rep2_trim0.01"
## [1] "omp_hmc_rep1_mean_omp_hmc_rep2_trim0.01_mean"
## [1] "d3xog_het_hmc_paired_q30_trim0.01"
## [1] "d3xog_het_hmc_paired_q30_trim0.01_mean"
## [1] "d3xog_ko_hmc_paired_q30_trim0.01"
## [1] "d3xog_ko_hmc_paired_q30_trim0.01_mean"
```

![plot of chunk gene_whole_W200N50F50_d3xog_hmc](figure/gene_whole_W200N50F50_d3xog_hmc.png) 

```
## [1] 0.0 0.8
```



```r
plot2.several("gene_whole_W200N50F50_chr", "d3xog_mc", data_type = "rpkm/mean", 
    cols = col3, group2 = "trim0.01")
```

```
## [1] "omp_mc_rmdup_trim0.01"
## [1] "omp_mc_rmdup_trim0.01_mean"
## [1] "d3xog_het_mc_paired_q30_trim0.01"
## [1] "d3xog_het_mc_paired_q30_trim0.01_mean"
## [1] "d3xog_ko_mc_paired_q30_trim0.01"
## [1] "d3xog_ko_mc_paired_q30_trim0.01_mean"
```

![plot of chunk gene_whole_W200N50F50_d3xog_mc](figure/gene_whole_W200N50F50_d3xog_mc.png) 

```
## [1] 0.080 0.759
```





```r
positionMatrix.all("gene_whole_W200N50F50_chr", data_type = "rpkm/mean")
```



```r
wt.hmc <- makeImage("omp_hmc_rep1_mean_omp_hmc_rep2", "gene_whole_W200N50F50_chr", 
    data_type = "rpkm/mean", image = F)
```

```
## [1] "/media/storage2/analysis/profiles/norm/rpkm/mean/gene_whole_W200N50F50_chr/images/omp_hmc_rep1_mean_omp_hmc_rep2"
```

```r
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



```r
wt.hmc.pc <- prcomp(wt.hmc[, 51:100])
wt.hmc.pred <- predict(wt.hmc.pc, wt.hmc)
```



```r
MP.heat(wt.hmc[order(wt.hmc.pred[, 1]), ], range = c(0, 2), average = 50)
```

![plot of chunk unnamed-chunk-10](figure/unnamed-chunk-10.png) 



```r
MP.heat(het.hmc[order(wt.hmc.pred[, 1]), ], range = c(0, 2), average = 50)
```

![plot of chunk unnamed-chunk-11](figure/unnamed-chunk-11.png) 



```r
MP.heat(ko.hmc[order(wt.hmc.pred[, 1]), ], range = c(0, 2), average = 50)
```

![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-12.png) 



