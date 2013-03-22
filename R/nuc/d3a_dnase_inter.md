D3a - DNase peaks
========================================================


```r
suppressPackageStartupMessages(source("~/src/seqAnalysis/R/profiles2.R"))
```


#### 5hmC

```r
plot2.several("d3a_het_dnase_sort_q30_dnase_sub_igenome_ensembl_genes_nosex.bed_W25F200_both_chr", 
    "d3a_hmc", data_type = "rpkm/mean", group2 = "trim0.01", cols = col2, y.vals = c(0, 
        1))
```

```
## [1] "moe_d3a_wt_hmc_rpkm_trim0.01"
## [1] "moe_d3a_wt_hmc_rpkm_trim0.01_mean"
## [1] "moe_d3a_ko_hmc_rpkm_trim0.01"
## [1] "moe_d3a_ko_hmc_rpkm_trim0.01_mean"
```

![plot of chunk unnamed-chunk-2](figure/unnamed-chunk-2.png) 

```
## [1] 0 1
```


#### 5mC

```r
plot2.several("d3a_het_dnase_sort_q30_dnase_sub_igenome_ensembl_genes_nosex.bed_W25F200_both_chr", 
    "d3a_mc", data_type = "rpkm/mean", group2 = "trim0.01", cols = col2, y.vals = c(0, 
        0.6))
```

```
## [1] "moe_d3a_wt_mc_rpkm_trim0.01"
## [1] "moe_d3a_wt_mc_rpkm_trim0.01_mean"
## [1] "moe_d3a_ko_mc_rpkm_trim0.01"
## [1] "moe_d3a_ko_mc_rpkm_trim0.01_mean"
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3.png) 

```
## [1] 0.0 0.6
```


#### Nucleosomes

```r
plot2.several("d3a_het_dnase_sort_q30_dnase_sub_igenome_ensembl_genes_nosex.bed_W25F200_both_chr", 
    "d3a_nuc_extend", data_type = "bam_extend/mean_chrom_mean", group2 = "trim0.01", 
    cols = col2, y.vals = c(0.5, 1.1))
```

```
## [1] "d3xog_wt_nuc_478_rmdup_q30_extend_trim0.01"
## [1] "d3xog_wt_nuc_478_rmdup_q30_extend_trim0.01_mean"
## [1] "d3xog_ko_nuc_256_rmdup_q30_extend_trim0.01"
## [1] "d3xog_ko_nuc_256_rmdup_q30_extend_trim0.01_mean"
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4.png) 

```
## [1] 0.5 1.1
```


#### DNase

```r
plot2("d3a_het_dnase_sort_q30_dnase_sub_igenome_ensembl_genes_nosex.bed_W25F200_both_chr", 
    "d3a_het_dnase_sort_q30", data_type = "bam_ends/mean", group2 = "trim0.01", 
    y.vals = c(0, 20))
```

```
## [1] "d3a_het_dnase_sort_q30_trim0.01"
## [1] "d3a_het_dnase_sort_q30_trim0.01_mean"
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5.png) 

```
## [1]  0 20
```


