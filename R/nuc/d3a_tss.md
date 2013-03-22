D3a TSSs
========================================================

```r
suppressPackageStartupMessages(source("~/src/seqAnalysis/R/profiles2.R"))
```


#### 5hmC

```r
plot2.several("refGene_noRandom_order_outsides2_tss_W25F200_chr", "d3a_hmc", 
    data_type = "rpkm/mean", group2 = "exclude-chrX_trim0.01", cols = col2)
```

```
## [1] "moe_d3a_wt_hmc_rpkm_exclude-chrX_trim0.01"
## [1] "moe_d3a_wt_hmc_rpkm_exclude-chrX_trim0.01_mean"
## [1] "moe_d3a_ko_hmc_rpkm_exclude-chrX_trim0.01"
## [1] "moe_d3a_ko_hmc_rpkm_exclude-chrX_trim0.01_mean"
```

![plot of chunk unnamed-chunk-2](figure/unnamed-chunk-2.png) 

```
## [1] 0.043 0.842
```


#### 5mC

```r
plot2.several("refGene_noRandom_order_outsides2_tss_W25F200_chr", "d3a_mc", 
    data_type = "rpkm/mean", group2 = "exclude-chrX_trim0.01", cols = col2)
```

```
## [1] "moe_d3a_wt_mc_rpkm_exclude-chrX_trim0.01"
## [1] "moe_d3a_wt_mc_rpkm_exclude-chrX_trim0.01_mean"
## [1] "moe_d3a_ko_mc_rpkm_exclude-chrX_trim0.01"
## [1] "moe_d3a_ko_mc_rpkm_exclude-chrX_trim0.01_mean"
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3.png) 

```
## [1] 0.057 0.553
```


#### Nucleosomes

```r
plot2.several("refGene_noRandom_order_outsides2_tss_W25F200_chr", "d3a_nuc_extend", 
    data_type = "bam_extend/mean_chrom_mean", group2 = "exclude-chrX_trim0.01", 
    cols = col2)
```

```
## [1] "d3xog_wt_nuc_478_rmdup_q30_extend_exclude-chrX_trim0.01"
## [1] "d3xog_wt_nuc_478_rmdup_q30_extend_exclude-chrX_trim0.01_mean"
## [1] "d3xog_ko_nuc_256_rmdup_q30_extend_exclude-chrX_trim0.01"
## [1] "d3xog_ko_nuc_256_rmdup_q30_extend_exclude-chrX_trim0.01_mean"
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4.png) 

```
## [1] 0.554 1.040
```



```r
makeProfile2.allSamp("refGene_noRandom_order_outsides2_tss_W25F200_chr", data_type = "bam_extend/mean_chrom_mean", 
    group2 = "d3xog_wt_rmrna_fpkm_range4", rm.outlier = 0.01)
```

```
## [1] "/media/storage2/analysis/profiles/norm/bam_extend/mean_chrom_mean/refGene_noRandom_order_outsides2_tss_W25F200_chr"
## Note: next may be used in wrong context: no loop is visible 
## [1] "Skipping"
## [1] "Skipping"
```

```
## Error: task 1 failed - "no loop for break/next, jumping to top level"
```

```r
makeProfile2.allSamp("refGene_noRandom_order_outsides2_tss_W25F200_chr", data_type = "bam_extend/mean_chrom_mean", 
    group2 = "d3xog_ko_rmrna_fpkm_range4", rm.outlier = 0.01)
```

```
## [1] "/media/storage2/analysis/profiles/norm/bam_extend/mean_chrom_mean/refGene_noRandom_order_outsides2_tss_W25F200_chr"
## Note: next may be used in wrong context: no loop is visible 
## [1] "Skipping"
## [1] "Skipping"
```

```
## Error: task 1 failed - "no loop for break/next, jumping to top level"
```



```r
plot2("refGene_noRandom_order_outsides2_tss_W25F200_chr", "d3xog_wt_nuc_478_rmdup_q30_extend", 
    group2 = "d3xog_wt_rmrna_fpkm_range4_trim0.01", data_type = "bam_extend/mean_chrom_mean", 
    cols = col4, y.vals = c(0.2, 1.2))
```

```
## [1] "d3xog_wt_nuc_478_rmdup_q30_extend_d3xog_wt_rmrna_fpkm_range4_trim0.01"
## [1] "d3xog_wt_nuc_478_rmdup_q30_extend_d3xog_wt_rmrna_fpkm_range4_trim0.01_mean"
```

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6.png) 

```
## [1] 0.2 1.2
```



```r
plot2("refGene_noRandom_order_outsides2_tss_W25F200_chr", "d3xog_ko_nuc_256_rmdup_q30_extend", 
    group2 = "d3xog_ko_rmrna_fpkm_range4_trim0.01", data_type = "bam_extend/mean_chrom_mean", 
    cols = col4, y.vals = c(0.2, 1.2))
```

```
## [1] "d3xog_ko_nuc_256_rmdup_q30_extend_d3xog_ko_rmrna_fpkm_range4_trim0.01"
## [1] "d3xog_ko_nuc_256_rmdup_q30_extend_d3xog_ko_rmrna_fpkm_range4_trim0.01_mean"
```

![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-7.png) 

```
## [1] 0.2 1.2
```



```r
par(mfrow = c(2, 2), mar = c(2, 2, 2, 2))
for (i in 1:4) {
    plot2.several("refGene_noRandom_order_outsides2_tss_W25F200_chr", "d3a_nuc_extend", 
        data_type = "bam_extend/mean_chrom_mean", group2 = "d3xog_wt_rmrna_fpkm_range4_trim0.01", 
        cols = col2, group2_col = i, y.vals = c(0.4, 1.2), lab = c("TSS"))
}
```

```
## [1] "d3xog_wt_nuc_478_rmdup_q30_extend_d3xog_wt_rmrna_fpkm_range4_trim0.01"
## [1] "d3xog_wt_nuc_478_rmdup_q30_extend_d3xog_wt_rmrna_fpkm_range4_trim0.01_mean"
## [1] "d3xog_ko_nuc_256_rmdup_q30_extend_d3xog_wt_rmrna_fpkm_range4_trim0.01"
## [1] "d3xog_ko_nuc_256_rmdup_q30_extend_d3xog_wt_rmrna_fpkm_range4_trim0.01_mean"
```

```
## [1] 0.4 1.2
## [1] "d3xog_wt_nuc_478_rmdup_q30_extend_d3xog_wt_rmrna_fpkm_range4_trim0.01"
## [1] "d3xog_wt_nuc_478_rmdup_q30_extend_d3xog_wt_rmrna_fpkm_range4_trim0.01_mean"
## [1] "d3xog_ko_nuc_256_rmdup_q30_extend_d3xog_wt_rmrna_fpkm_range4_trim0.01"
## [1] "d3xog_ko_nuc_256_rmdup_q30_extend_d3xog_wt_rmrna_fpkm_range4_trim0.01_mean"
```

```
## [1] 0.4 1.2
## [1] "d3xog_wt_nuc_478_rmdup_q30_extend_d3xog_wt_rmrna_fpkm_range4_trim0.01"
## [1] "d3xog_wt_nuc_478_rmdup_q30_extend_d3xog_wt_rmrna_fpkm_range4_trim0.01_mean"
## [1] "d3xog_ko_nuc_256_rmdup_q30_extend_d3xog_wt_rmrna_fpkm_range4_trim0.01"
## [1] "d3xog_ko_nuc_256_rmdup_q30_extend_d3xog_wt_rmrna_fpkm_range4_trim0.01_mean"
```

```
## [1] 0.4 1.2
## [1] "d3xog_wt_nuc_478_rmdup_q30_extend_d3xog_wt_rmrna_fpkm_range4_trim0.01"
## [1] "d3xog_wt_nuc_478_rmdup_q30_extend_d3xog_wt_rmrna_fpkm_range4_trim0.01_mean"
## [1] "d3xog_ko_nuc_256_rmdup_q30_extend_d3xog_wt_rmrna_fpkm_range4_trim0.01"
## [1] "d3xog_ko_nuc_256_rmdup_q30_extend_d3xog_wt_rmrna_fpkm_range4_trim0.01_mean"
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8.png) 

```
## [1] 0.4 1.2
```


#### DNase

```r
plot2("refGene_noRandom_order_outsides2_tss_W25F200_chr", "d3a_het_dnase_sort_q30", 
    data_type = "bam_ends/mean", group2 = "exclude-chrX_trim0.01")
```

```
## [1] "d3a_het_dnase_sort_q30_exclude-chrX_trim0.01"
## [1] "d3a_het_dnase_sort_q30_exclude-chrX_trim0.01_mean"
```

![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-9.png) 

```
## [1] 1.398 7.795
```

