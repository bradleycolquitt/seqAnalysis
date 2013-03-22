D3a - DNase peaks 2 
========================================================


```r
suppressPackageStartupMessages(source("~/src/seqAnalysis/R/profiles2.R"))
```



```r
makeProfile2.allSamp("d3a_het_dnase_sort_q30_dnase_sub_igenome_ensembl_genes_extend5kb.bed_W25F200_both_chr", 
    data_type = "bam_extend/mean_chrom_mean", rm.outliers = 0.01)
makeProfile2.allSamp("d3a_het_dnase_sort_q30_dnase_sub_igenome_ensembl_genes_extend5kb.bed_W25F200_both_chr", 
    data_type = "bam_ends/mean", rm.outliers = 0.01)
makeProfile2.allSamp("d3a_het_dnase_sort_q30_dnase_sub_igenome_ensembl_genes_extend5kb.bed_W25F200_both_chr", 
    data_type = "rpkm/mean", rm.outliers = 0.01)
```



```r
makeProfile2.allSamp("d3a_het_dnase_sort_q30_dnase_sub_igenome_ensembl_genes_extend5kb_shuffle.bed_W25F200_both_chr", 
    data_type = "bam_extend/mean_chrom_mean", rm.outliers = 0.01)
makeProfile2.allSamp("d3a_het_dnase_sort_q30_dnase_sub_igenome_ensembl_genes_extend5kb_shuffle.bed_W25F200_both_chr", 
    data_type = "bam_ends/mean", rm.outliers = 0.01)
makeProfile2.allSamp("d3a_het_dnase_sort_q30_dnase_sub_igenome_ensembl_genes_extend5kb_shuffle.bed_W25F200_both_chr", 
    data_type = "rpkm/mean", rm.outliers = 0.01)
```


#### 5hmC

```r
par(mfrow = c(1, 2), mar = c(2, 2, 2, 2))
plot2.several("d3a_het_dnase_sort_q30_dnase_sub_igenome_ensembl_genes_extend5kb.bed_W25F200_both_chr", 
    "d3a_hmc", data_type = "rpkm/mean", group2 = "trim0.01", cols = col2, y.vals = c(0, 
        1))
```

```
## [1] "moe_d3a_wt_hmc_rpkm_trim0.01"
## [1] "moe_d3a_wt_hmc_rpkm_trim0.01_mean"
## [1] "moe_d3a_ko_hmc_rpkm_trim0.01"
## [1] "moe_d3a_ko_hmc_rpkm_trim0.01_mean"
```

```
## [1] 0 1
```

```r
plot2.several("d3a_het_dnase_sort_q30_dnase_sub_igenome_ensembl_genes_extend5kb_shuffle.bed_W25F200_both_chr", 
    "d3a_hmc", data_type = "rpkm/mean", group2 = "trim0.01", cols = col2, y.vals = c(0, 
        1))
```

```
## [1] "moe_d3a_wt_hmc_rpkm_trim0.01"
## [1] "moe_d3a_wt_hmc_rpkm_trim0.01_mean"
## [1] "moe_d3a_ko_hmc_rpkm_trim0.01"
## [1] "moe_d3a_ko_hmc_rpkm_trim0.01_mean"
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4.png) 

```
## [1] 0 1
```


#### 5mC

```r
par(mfrow = c(1, 2), mar = c(2, 2, 2, 2))
plot2.several("d3a_het_dnase_sort_q30_dnase_sub_igenome_ensembl_genes_extend5kb.bed_W25F200_both_chr", 
    "d3a_mc", data_type = "rpkm/mean", group2 = "trim0.01", cols = col2, y.vals = c(0, 
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
plot2.several("d3a_het_dnase_sort_q30_dnase_sub_igenome_ensembl_genes_extend5kb_shuffle.bed_W25F200_both_chr", 
    "d3a_hmc", data_type = "rpkm/mean", group2 = "trim0.01", cols = col2, y.vals = c(0, 
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
## [1] 0 1
```


#### Nucleosomes

```r
plot2.several("d3a_het_dnase_sort_q30_dnase_sub_igenome_ensembl_genes_extend5kb.bed_W25F200_both_chr", 
    "d3a_nuc_extend", data_type = "bam_extend/mean_chrom_mean", group2 = "trim0.01", 
    cols = col2, y.vals = c(0.7, 1.1))
```

```
## [1] "d3xog_wt_nuc_478_rmdup_q30_extend_trim0.01"
## [1] "d3xog_wt_nuc_478_rmdup_q30_extend_trim0.01_mean"
## [1] "d3xog_ko_nuc_256_rmdup_q30_extend_trim0.01"
## [1] "d3xog_ko_nuc_256_rmdup_q30_extend_trim0.01_mean"
```

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6.png) 

```
## [1] 0.7 1.1
```


Random genomic sites

```r
plot2.several("d3a_het_dnase_sort_q30_dnase_sub_igenome_ensembl_genes_extend5kb_shuffle.bed_W25F200_both_chr", 
    "d3a_nuc_extend", data_type = "bam_extend/mean_chrom_mean", group2 = "trim0.01", 
    cols = col2, y.vals = c(0.7, 1.1))
```

```
## [1] "d3xog_wt_nuc_478_rmdup_q30_extend_trim0.01"
## [1] "d3xog_wt_nuc_478_rmdup_q30_extend_trim0.01_mean"
## [1] "d3xog_ko_nuc_256_rmdup_q30_extend_trim0.01"
## [1] "d3xog_ko_nuc_256_rmdup_q30_extend_trim0.01_mean"
```

![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-7.png) 

```
## [1] 0.7 1.1
```



#### DNase

```r
plot2("d3a_het_dnase_sort_q30_dnase_sub_igenome_ensembl_genes_extend5kb.bed_W25F200_both_chr", 
    "d3a_het_dnase_sort_q30", data_type = "bam_ends/mean", group2 = "trim0.01", 
    y.vals = c(0, 0.7))
```

```
## [1] "d3a_het_dnase_sort_q30_trim0.01"
## [1] "d3a_het_dnase_sort_q30_trim0.01_mean"
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8.png) 

```
## [1] 0.0 0.7
```


Expression of genes adjacent to intergenic DNase HS sites with increased nuc occupancy
------------------

```r
d3.inter.q75.rg <- read.delim("~/s2/data/homer/peaks/d3a_het_dnase_sort_q30/d3a_het_dnase_sort_q30_dnase_sub_igenome_ensembl_genes_extend5kb_nosex_mid1kb_d3xog_ko_wt_ratio_q75_closest_refgene_nodup.bed", 
    header = F)
rna.1log2 <- readRDS("~/s2/analysis/rna/rdata/d3xog_wt_ko_rmrna_masked_uq_comp_1log2.rds")
rna.1log2$d3.inter.q75 <- FALSE
rna.1log2$d3.inter.q75[rna.1log2$gene %in% d3.inter.q75.rg[, 10]] <- TRUE
```



```r
theme_set(theme_bw())
```

```
## Error: could not find function "theme_set"
```

```r
gg <- ggplot(rna.1log2, aes(x = ko.wt, color = d3.inter.q75))
```

```
## Error: could not find function "ggplot"
```

```r
gg <- gg + geom_density() + coord_cartesian(xlim = c(-3, 3)) + scale_color_manual(values = c("black", 
    "red"))
```

```
## Error: object 'gg' not found
```

```r
gg <- gg + labs(x = c("KO / WT log2(FPKM + 1)"))
```

```
## Error: object 'gg' not found
```

```r
gg
```

```
## Error: object 'gg' not found
```

