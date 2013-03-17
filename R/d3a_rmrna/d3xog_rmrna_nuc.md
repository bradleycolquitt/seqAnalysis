D3xOG rmRNA and nucleosomes
========================================================


```r
suppressPackageStartupMessages(source("~/src/seqAnalysis/R/image.R"))
suppressPackageStartupMessages(source("~/src/seqAnalysis/R/profiles2.R"))
nuc.wt <- makeImage("d3xog_wt_nuc_478", "refGene_noRandom_order_outsides2_tss_W25F200_chr", 
    data_type = "rpkm/mean", image = F)
```

```
## [1] "/media/storage2/analysis/profiles/norm/rpkm/mean/refGene_noRandom_order_outsides2_tss_W25F200_chr/images/d3xog_wt_nuc_478"
```

```r
nuc.ko <- makeImage("d3xog_ko_nuc_256", "refGene_noRandom_order_outsides2_tss_W25F200_chr", 
    data_type = "rpkm/mean", image = F)
```

```
## [1] "/media/storage2/analysis/profiles/norm/rpkm/mean/refGene_noRandom_order_outsides2_tss_W25F200_chr/images/d3xog_ko_nuc_256"
```

```r

nuc.na <- which(is.na(nuc.wt), arr.ind = TRUE)
nuc.wt[nuc.na] <- 0

nuc.na <- which(is.na(nuc.ko), arr.ind = TRUE)
nuc.ko[nuc.na] <- 0

nuc.ko.wt <- nuc.ko - nuc.wt
```



```r
rna <- readRDS("~/s2/analysis/rna/rdata/d3xog_wt_ko_rmrna_masked_uq_comp_1log2.rds")

nuc.ko.wt.ord.rna <- nuc.ko.wt[match(rna$gene[order(rna$ko.wt)], rownames(nuc.ko.wt)), 
    ]
nuc.ko.wt.ord.rna <- na.omit(nuc.ko.wt.ord.rna)
```



```r
MP.heat(nuc.ko.wt, range = c(-1, 1), average = 50)
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3.png) 



```r
nuc.ko.wt.ord.rna.down <- apply(nuc.ko.wt.ord.rna[, 201:220], 1, mean)
nuc.ko.wt.ord.rna.down.c100 <- foreach(c = isplitVector(nuc.ko.wt.ord.rna.down, 
    chunks = 100), .combine = "c") %do% mean(c, na.rm = TRUE)
nuc.ko.wt.ord.rna.down.c100 <- as.data.frame(nuc.ko.wt.ord.rna.down.c100)
nuc.ko.wt.ord.rna.down.c100$index <- 100:1
nuc.ko.wt.ord.rna.down.c100.boot <- foreach(c = isplitVector(nuc.ko.wt.ord.rna.down, 
    chunks = 100), .combine = "rbind") %do% bootCI(c)
nuc.ko.wt.ord.rna.down.c100 <- cbind(nuc.ko.wt.ord.rna.down.c100, nuc.ko.wt.ord.rna.down.c100.boot)
colnames(nuc.ko.wt.ord.rna.down.c100)[3:4] <- c("lower", "upper")
```



```r
nuc.ko.wt.ord.rna.down.c100$wilcox.FDR <- p.adjust(foreach(c = isplitVector(nuc.ko.wt.ord.rna.down, 
    chunks = 100), .combine = "c") %do% wilcox.test(c)$p.value, method = "fdr")
nuc.ko.wt.ord.rna.down.c100$wilcox.FDR.05 <- cut(nuc.ko.wt.ord.rna.down.c100$wilcox.FDR, 
    breaks = c(0, 0.05, 1))
```



```r
theme_set(theme_gray())
```

```
## Error: could not find function "theme_set"
```

```r
gg <- ggplot(nuc.ko.wt.ord.rna.down.c100, aes(nuc.ko.wt.ord.rna.down.c100, index))
```

```
## Error: could not find function "ggplot"
```

```r
gg <- gg + geom_vline(xintercept = 0, color = "red")
```

```
## Error: object 'gg' not found
```

```r
gg <- gg + geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0, size = 0.1) + 
    geom_point(aes(color = wilcox.FDR.05), size = 2) + xlab("RPM") + ylab("") + 
    theme(legend.position = "none", axis.text.y = element_blank()) + labs(title = c("KO - WT 5hmC")) + 
    scale_color_manual(values = c("red", "black"))
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


**General reduction of nucleosome occupancy from TSS to 1kb downstream, irrspective of transcriptional change**

### Profiles

```r
suppressPackageStartupMessages(source("~/src/seqAnalysis/R/profiles2.R"))
```


#### Top expressors
**Nucleosomes**

```r
plot2.several("refGene_noRandom_order_outsides2_tss_W25F200_chr", "d3a_nuc2", 
    group2 = "d3xog_ko_rmrna_fpkm_range4", data_type = "rpkm/mean", cols = col2, 
    fname = "manual", group2_col = 4, lab = c("TSS"))
```

```
## [1] "d3xog_wt_nuc_478_rmdup_d3xog_ko_rmrna_fpkm_range4"
## [1] "d3xog_wt_nuc_478_rmdup_d3xog_ko_rmrna_fpkm_range4_mean"
## [1] "d3xog_ko_nuc_256_rmdup_d3xog_ko_rmrna_fpkm_range4"
## [1] "d3xog_ko_nuc_256_rmdup_d3xog_ko_rmrna_fpkm_range4_mean"
```

![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-7.png) 

```
## [1] 0.053 0.538
```


**5hmC**

```r
a <- plot2.several("refGene_noRandom_order_outsides2_tss_W25F200_chr", "d3a_hmc", 
    group2 = "moe_d3a_wt_mrna_range4", data_type = "rpkm/mean", cols = col2, 
    fname = "manual", group2_col = 4, lab = c("TSS"))
```

```
## [1] "moe_d3a_wt_hmc_rpkm_moe_d3a_wt_mrna_range4"
## [1] "moe_d3a_wt_hmc_rpkm_moe_d3a_wt_mrna_range4_mean"
## [1] "moe_d3a_ko_hmc_rpkm_moe_d3a_wt_mrna_range4"
## [1] "moe_d3a_ko_hmc_rpkm_moe_d3a_wt_mrna_range4_mean"
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8.png) 

```
## [1] -0.202  1.398
```


**5mC**

```r
a <- plot2.several("refGene_noRandom_order_outsides2_tss_W25F200_chr", "d3a_mc", 
    group2 = "moe_d3a_wt_mrna_range4", data_type = "rpkm/mean", cols = col2, 
    fname = "manual", group2_col = 4, lab = c("TSS"))
```

```
## [1] "moe_d3a_wt_mc_rpkm_moe_d3a_wt_mrna_range4"
## [1] "moe_d3a_wt_mc_rpkm_moe_d3a_wt_mrna_range4_mean"
## [1] "moe_d3a_ko_mc_rpkm_moe_d3a_wt_mrna_range4"
## [1] "moe_d3a_ko_mc_rpkm_moe_d3a_wt_mrna_range4_mean"
```

![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-9.png) 

```
## [1] -0.073  0.659
```


#### Silent (FPKM = 0)
**Nucleosomes**

```r
plot2.several("refGene_noRandom_order_outsides2_tss_W25F200_chr", "d3a_nuc2", 
    group2 = "d3xog_ko_rmrna_fpkm_range4", data_type = "rpkm/mean", cols = col2, 
    fname = "manual", group2_col = 1, lab = c("TSS"))
```

```
## [1] "d3xog_wt_nuc_478_rmdup_d3xog_ko_rmrna_fpkm_range4"
## [1] "d3xog_wt_nuc_478_rmdup_d3xog_ko_rmrna_fpkm_range4_mean"
## [1] "d3xog_ko_nuc_256_rmdup_d3xog_ko_rmrna_fpkm_range4"
## [1] "d3xog_ko_nuc_256_rmdup_d3xog_ko_rmrna_fpkm_range4_mean"
```

![plot of chunk unnamed-chunk-10](figure/unnamed-chunk-10.png) 

```
## [1] 0.297 0.611
```


**5hmC**

```r
a <- plot2.several("refGene_noRandom_order_outsides2_tss_W25F200_chr", "d3a_hmc", 
    group2 = "moe_d3a_wt_mrna_range4", data_type = "rpkm/mean", cols = col2, 
    fname = "manual", group2_col = 1, lab = c("TSS"))
```

```
## [1] "moe_d3a_wt_hmc_rpkm_moe_d3a_wt_mrna_range4"
## [1] "moe_d3a_wt_hmc_rpkm_moe_d3a_wt_mrna_range4_mean"
## [1] "moe_d3a_ko_hmc_rpkm_moe_d3a_wt_mrna_range4"
## [1] "moe_d3a_ko_hmc_rpkm_moe_d3a_wt_mrna_range4_mean"
```

![plot of chunk unnamed-chunk-11](figure/unnamed-chunk-11.png) 

```
## [1] 0.191 0.495
```


**5mC**

```r
a <- plot2.several("refGene_noRandom_order_outsides2_tss_W25F200_chr", "d3a_mc", 
    group2 = "moe_d3a_wt_mrna_range4", data_type = "rpkm/mean", cols = col2, 
    fname = "manual", group2_col = 1, lab = c("TSS"))
```

```
## [1] "moe_d3a_wt_mc_rpkm_moe_d3a_wt_mrna_range4"
## [1] "moe_d3a_wt_mc_rpkm_moe_d3a_wt_mrna_range4_mean"
## [1] "moe_d3a_ko_mc_rpkm_moe_d3a_wt_mrna_range4"
## [1] "moe_d3a_ko_mc_rpkm_moe_d3a_wt_mrna_range4_mean"
```

![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-12.png) 

```
## [1] 0.181 0.682
```

