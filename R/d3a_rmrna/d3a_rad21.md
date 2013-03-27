D3a Rad21 sites
========================================================


```r
suppressPackageStartupMessages(source("~/src/seqAnalysis/R/profiles2.R"))
```




```r
makeProfile2.allSamp("cad_n2a_rad21_r12_rand10E4_nosex.bed_W25F200_both_chr", 
    data_type = "bam_extend/mean_chrom_mean_0", rm.outliers = 0.01)
makeProfile2.allSamp("cad_n2a_rad21_r12_rand10E4_nosex.bed_W25F200_both_chr", 
    data_type = "bam_ends/mean", rm.outliers = 0.01)
makeProfile2.allSamp("cad_n2a_rad21_r12_rand10E4_nosex.bed_W25F200_both_chr", 
    data_type = "rpkm/mean", rm.outliers = 0.01)
```



```r
makeProfile2.allSamp("cad_n2a_rad21_r12_rand10E4_nosex.bed_W200F500_both_chr", 
    data_type = "bam_extend/mean_chrom_mean_0", rm.outliers = 0.01)
makeProfile2.allSamp("cad_n2a_rad21_r12_rand10E4_nosex.bed_W200F500_both_chr", 
    data_type = "bam_ends/mean", rm.outliers = 0.01)
makeProfile2.allSamp("cad_n2a_rad21_r12_rand10E4_nosex.bed_W200F500_both_chr", 
    data_type = "rpkm/mean", rm.outliers = 0.01)
```



```r
makeProfile2.allSamp("cad_n2a_rad21_r12_rand10E4_shuffle_nosex.bed_W25F200_both_chr", 
    data_type = "bam_extend/mean_chrom_mean_0", rm.outliers = 0.01)
makeProfile2.allSamp("cad_n2a_rad21_r12_rand10E4_shuffle_nosex.bed_W25F200_both_chr", 
    data_type = "bam_ends/mean", rm.outliers = 0.01)
makeProfile2.allSamp("cad_n2a_rad21_r12_rand10E4_shuffle_nosex.bed_W25F200_both_chr", 
    data_type = "rpkm/mean", rm.outliers = 0.01)
```


#### 5hmC

```r
par(mfrow = c(1, 2), mar = c(2, 2, 2, 2))
plot2.several("cad_n2a_rad21_r12_rand10E4_nosex.bed_W25F200_both_chr", "d3a_hmc", 
    data_type = "rpkm/mean", group2 = "trim0.01", cols = col2, y.vals = c(0.2, 
        1))
```

```
## [1] "moe_d3a_wt_hmc_rpkm_trim0.01"
## [1] "moe_d3a_wt_hmc_rpkm_trim0.01_mean"
## [1] "moe_d3a_ko_hmc_rpkm_trim0.01"
## [1] "moe_d3a_ko_hmc_rpkm_trim0.01_mean"
```

```
## [1] 0.2 1.0
```

```r
plot2.several("cad_n2a_rad21_r12_rand10E4_shuffle_nosex.bed_W25F200_both_chr", 
    "d3a_hmc", data_type = "rpkm/mean", group2 = "trim0.01", cols = col2, y.vals = c(0.2, 
        1))
```

```
## [1] "moe_d3a_wt_hmc_rpkm_trim0.01"
## [1] "moe_d3a_wt_hmc_rpkm_trim0.01_mean"
## [1] "moe_d3a_ko_hmc_rpkm_trim0.01"
## [1] "moe_d3a_ko_hmc_rpkm_trim0.01_mean"
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5.png) 

```
## [1] 0.2 1.0
```



```r
par(mfrow = c(1, 1), mar = c(2, 2, 2, 2))
plot2.several("cad_n2a_rad21_r12_rand10E4_nosex.bed_W200F500_both_chr", "d3a_hmc", 
    data_type = "rpkm/mean", group2 = "trim0.01", cols = col2, y.vals = c(0.2, 
        1), range = c(251, 750))
```

```
## [1] "moe_d3a_wt_hmc_rpkm_trim0.01"
## [1] "moe_d3a_wt_hmc_rpkm_trim0.01_mean"
## [1] "moe_d3a_ko_hmc_rpkm_trim0.01"
## [1] "moe_d3a_ko_hmc_rpkm_trim0.01_mean"
```

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6.png) 

```
## [1] 0.2 1.0
```

```r
#
# plot2.several('cad_n2a_rad21_r12_rand10E4_shuffle_nosex.bed_W25F200_both_chr',
# 'd3a_hmc', data_type='rpkm/mean', group2='trim0.01', cols=col2,
# y.vals=c(.2, 1))
```

#### 5mC

```r
par(mfrow = c(1, 2), mar = c(2, 2, 2, 2))
plot2.several("cad_n2a_rad21_r12_rand10E4_nosex.bed_W25F200_both_chr", "d3a_mc", 
    data_type = "rpkm/mean", group2 = "trim0.01", cols = col2, y.vals = c(0, 
        0.6))
```

```
## [1] "moe_d3a_wt_mc_rpkm_trim0.01"
## [1] "moe_d3a_wt_mc_rpkm_trim0.01_mean"
## [1] "moe_d3a_ko_mc_rpkm_trim0.01"
## [1] "moe_d3a_ko_mc_rpkm_trim0.01_mean"
```

```
## [1] 0.0 0.6
```

```r
plot2.several("cad_n2a_rad21_r12_rand10E4_shuffle_nosex.bed_W25F200_both_chr", 
    "d3a_mc", data_type = "rpkm/mean", group2 = "trim0.01", cols = col2, y.vals = c(0, 
        0.6))
```

```
## [1] "moe_d3a_wt_mc_rpkm_trim0.01"
## [1] "moe_d3a_wt_mc_rpkm_trim0.01_mean"
## [1] "moe_d3a_ko_mc_rpkm_trim0.01"
## [1] "moe_d3a_ko_mc_rpkm_trim0.01_mean"
```

![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-7.png) 

```
## [1] 0.0 0.6
```



```r
par(mfrow = c(1, 1), mar = c(2, 2, 2, 2))
plot2.several("cad_n2a_rad21_r12_rand10E4_nosex.bed_W200F500_both_chr", "d3a_mc", 
    data_type = "rpkm/mean", group2 = "trim0.01", cols = col2, y.vals = c(0, 
        0.6), range = c(251, 750))
```

```
## [1] "moe_d3a_wt_mc_rpkm_trim0.01"
## [1] "moe_d3a_wt_mc_rpkm_trim0.01_mean"
## [1] "moe_d3a_ko_mc_rpkm_trim0.01"
## [1] "moe_d3a_ko_mc_rpkm_trim0.01_mean"
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8.png) 

```
## [1] 0.0 0.6
```

```r
#
# plot2.several('cad_n2a_rad21_r12_rand10E4_shuffle_nosex.bed_W25F200_both_chr',
# 'd3a_mc', data_type='rpkm/mean', group2='trim0.01', cols=col2,
# y.vals=c(0, .6))
```


#### Nucleosomes

```r
par(mfrow = c(1, 2), mar = c(2, 2, 2, 2))
plot2.several("cad_n2a_rad21_r12_rand10E4_nosex.bed_W25F200_both_chr", "d3a_nuc_extend", 
    data_type = "bam_extend/mean_chrom_mean_0", group2 = "trim0.01", cols = col2, 
    y.vals = c(0.8, 1.3), baseline = F)
```

```
## [1] "d3xog_wt_nuc_478_rmdup_q30_extend_trim0.01"
## [1] "d3xog_wt_nuc_478_rmdup_q30_extend_trim0.01_mean"
## [1] "d3xog_ko_nuc_256_rmdup_q30_extend_trim0.01"
## [1] "d3xog_ko_nuc_256_rmdup_q30_extend_trim0.01_mean"
```

```
## [1] 0.8 1.3
```

```r
plot2.several("cad_n2a_rad21_r12_rand10E4_shuffle_nosex.bed_W25F200_both_chr", 
    "d3a_nuc_extend", data_type = "bam_extend/mean_chrom_mean_0", group2 = "trim0.01", 
    cols = col2, y.vals = c(0.8, 1.3), baseline = F)
```

```
## [1] "d3xog_wt_nuc_478_rmdup_q30_extend_trim0.01"
## [1] "d3xog_wt_nuc_478_rmdup_q30_extend_trim0.01_mean"
## [1] "d3xog_ko_nuc_256_rmdup_q30_extend_trim0.01"
## [1] "d3xog_ko_nuc_256_rmdup_q30_extend_trim0.01_mean"
```

![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-9.png) 

```
## [1] 0.8 1.3
```



```r
par(mfrow = c(1, 1), mar = c(2, 2, 2, 2))
plot2.several("cad_n2a_rad21_r12_rand10E4_nosex.bed_W200F500_both_chr", "d3a_nuc_extend", 
    data_type = "bam_extend/mean_chrom_mean_0", group2 = "trim0.01", cols = col2, 
    y.vals = c(0.8, 1.3), baseline = T, wsize = 200)
```

```
## [1] "d3xog_wt_nuc_478_rmdup_q30_extend_trim0.01"
## [1] "d3xog_wt_nuc_478_rmdup_q30_extend_trim0.01_mean"
## [1] "d3xog_ko_nuc_256_rmdup_q30_extend_trim0.01"
## [1] "d3xog_ko_nuc_256_rmdup_q30_extend_trim0.01_mean"
```

![plot of chunk unnamed-chunk-10](figure/unnamed-chunk-10.png) 

```
## [1] 0.8 1.3
```

```r
#
# plot2.several('cad_n2a_rad21_r12_rand10E4_shuffle_nosex.bed_W25F200_both_chr','d3a_nuc_extend',
# data_type='bam_extend/mean_chrom_mean_0', group2='trim0.01', cols=col2,
# y.vals=c(.8, 1.3), baseline=F)
```



```r
par(mfrow = c(1, 2))
plot2.several("cad_n2a_rad21_r12_rand10E4_nosex.bed_W25F200_both_chr", "d3a_nuc2", 
    data_type = "rpkm/mean", group2 = "trim0.01", cols = col2, baseline = T)
```

```
## [1] "d3xog_wt_nuc_478_rmdup_trim0.01"
## [1] "d3xog_wt_nuc_478_rmdup_trim0.01_mean"
## [1] "d3xog_ko_nuc_256_rmdup_trim0.01"
## [1] "d3xog_ko_nuc_256_rmdup_trim0.01_mean"
```

```
## [1] 0.296 0.473
```

```r

plot2.several("cad_n2a_rad21_r12_rand10E4_nosex.bed_W200F500_both_chr", "d3a_nuc2", 
    data_type = "rpkm/mean", group2 = "trim0.01", cols = col2, baseline = T, 
    wsize = 200, range = c(251, 750))
```

```
## [1] "d3xog_wt_nuc_478_rmdup_trim0.01"
## [1] "d3xog_wt_nuc_478_rmdup_trim0.01_mean"
## [1] "d3xog_ko_nuc_256_rmdup_trim0.01"
## [1] "d3xog_ko_nuc_256_rmdup_trim0.01_mean"
```

![plot of chunk unnamed-chunk-11](figure/unnamed-chunk-11.png) 

```
## [1] 0.350 0.466
```


#### DNase

```r
plot2("cad_n2a_rad21_r12_rand10E4_nosex.bed_W200F500_both_chr", "d3a_het_dnase_sort_q30", 
    data_type = "bam_ends/mean", group2 = "trim0.01", wsize = 200)
```

```
## [1] "d3a_het_dnase_sort_q30_trim0.01"
## [1] "d3a_het_dnase_sort_q30_trim0.01_mean"
```

![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-12.png) 

```
## [1] 0.370 0.797
```

