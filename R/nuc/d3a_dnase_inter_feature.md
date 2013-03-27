DNase HS sites - Dnmt3a 5hmC/5mC/Nuc
========================================================


```r
suppressPackageStartupMessages(source("~/src/seqAnalysis/R/features.R"))
library(plyr)
library(ggplot2)
```



```r
d0.mod0 <- makeFeatureMatrix2("d3a_het_dnase_sort_q30_dnase_sub_igenome_ensembl_genes_extend5kb_nosex.bed_chr", 
    "d3a_rpkm", data_type = "rpkm/mean")
```

```
## [1] "moe_d3a_wt_hmc_rpkm" "moe_d3a_ko_hmc_rpkm" "moe_d3a_wt_mc_rpkm" 
## [4] "moe_d3a_ko_mc_rpkm"
```

```r
d0.mod500 <- makeFeatureMatrix2("d3a_het_dnase_sort_q30_dnase_sub_igenome_ensembl_genes_extend5kb_nosex.bed_chr", 
    "d3a_rpkm", data_type = "rpkm/mean", select = "flank500")
```

```
##            moe_d3a_wt_hmc_rpkm            moe_d3a_ko_hmc_rpkm 
## "moe_d3a_wt_hmc_rpkm_flank500" "moe_d3a_ko_hmc_rpkm_flank500" 
##             moe_d3a_wt_mc_rpkm             moe_d3a_ko_mc_rpkm 
##  "moe_d3a_wt_mc_rpkm_flank500"  "moe_d3a_ko_mc_rpkm_flank500"
```

```r
d0.mod1000 <- makeFeatureMatrix2("d3a_het_dnase_sort_q30_dnase_sub_igenome_ensembl_genes_extend5kb_nosex.bed_chr", 
    "d3a_rpkm", data_type = "rpkm/mean", select = "flank1000")
```

```
##             moe_d3a_wt_hmc_rpkm             moe_d3a_ko_hmc_rpkm 
## "moe_d3a_wt_hmc_rpkm_flank1000" "moe_d3a_ko_hmc_rpkm_flank1000" 
##              moe_d3a_wt_mc_rpkm              moe_d3a_ko_mc_rpkm 
##  "moe_d3a_wt_mc_rpkm_flank1000"  "moe_d3a_ko_mc_rpkm_flank1000"
```

```r


d0.nuc0 <- makeFeatureMatrix2("d3a_het_dnase_sort_q30_dnase_sub_igenome_ensembl_genes_extend5kb_nosex.bed_chr", 
    "d3a_nuc_extend", data_type = "bam_extend/mean_chrom_mean_0")
```

```
## [1] "d3xog_wt_nuc_478_rmdup_q30_extend" "d3xog_ko_nuc_256_rmdup_q30_extend"
```

```r
d0.nuc500 <- makeFeatureMatrix2("d3a_het_dnase_sort_q30_dnase_sub_igenome_ensembl_genes_extend5kb_nosex.bed_chr", 
    "d3a_nuc_extend", data_type = "bam_extend/mean_chrom_mean_0", select = "flank500")
```

```
##            d3xog_wt_nuc_478_rmdup_q30_extend 
## "d3xog_wt_nuc_478_rmdup_q30_extend_flank500" 
##            d3xog_ko_nuc_256_rmdup_q30_extend 
## "d3xog_ko_nuc_256_rmdup_q30_extend_flank500"
```

```r
d0.nuc1000 <- makeFeatureMatrix2("d3a_het_dnase_sort_q30_dnase_sub_igenome_ensembl_genes_extend5kb_nosex.bed_chr", 
    "d3a_nuc_extend", data_type = "bam_extend/mean_chrom_mean_0", select = "flank1000")
```

```
##             d3xog_wt_nuc_478_rmdup_q30_extend 
## "d3xog_wt_nuc_478_rmdup_q30_extend_flank1000" 
##             d3xog_ko_nuc_256_rmdup_q30_extend 
## "d3xog_ko_nuc_256_rmdup_q30_extend_flank1000"
```

```r

```



```r
m <- match(rownames(d0.mod0), rownames(d0.nuc0))
d0.comb <- cbind(d0.mod0, d0.mod500, d0.mod1000, d0.nuc0[m, ], d0.nuc500[m, 
    ], d0.nuc1000[m, ])
d0.comb.melt <- melt(d0.comb)
```

```
## Using as id variables
```

```r
d0.comb.melt$id <- rownames(d0.comb)
d0.comb.melt$mod <- as.factor(c(rep(rep(c("hmc", "mc"), each = nrow(d0.mod0) * 
    2), times = 3), rep("nuc", each = nrow(d0.mod0) * 6)))
d0.comb.melt$genotype <- as.factor(rep(c("wt", "ko"), each = nrow(d0.mod0)))
d0.comb.melt$pos <- as.factor(c(rep(c("0", "500", "1000"), each = nrow(d0.mod0) * 
    4), rep(c("0", "500", "1000"), each = nrow(d0.mod0) * 2)))
```



```r
d0.comb.ratio <- ddply(d0.comb.melt, .(id, mod, pos), summarize, ko.wt = log2((value[genotype == 
    "ko"] + 0.01)/(value[genotype == "wt"] + 0.01)))
```



```r
d0.comb.ratio.c <- recast(d0.comb.ratio, id ~ mod + pos)
```

```
## Using id, mod, pos as id variables
```



```r
gg <- ggplot(d0.comb.ratio.c, aes(hmc_500, nuc_500))
gg + geom_point(alpha = I(1/5)) + labs(title = "5hmC")
```

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6.png) 



```r
gg <- ggplot(d0.comb.ratio.c, aes(mc_500, nuc_500))
gg + geom_point(alpha = I(1/5)) + labs(title = "5mC")
```

![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-7.png) 



```r
d0.comb.ratio.c$hmc_500_cut4 <- cut(d0.comb.ratio.c$hmc_500, quantile(d0.comb.ratio.c$hmc_500))
d0.comb.ratio$hmc_500_cut4 <- d0.comb.ratio.c$hmc_500_cut4[match(d0.comb.ratio$id, 
    d0.comb.ratio.c[, 1])]
```



```r
gg <- ggplot(na.omit(d0.comb.ratio[d0.comb.ratio$mod == "nuc" & d0.comb.ratio$pos == 
    "500", ]), aes(hmc_500_cut4, ko.wt))
gg + geom_boxplot() + coord_cartesian(ylim = c(-1, 2))
```

![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-9.png) 



```r
d0.comb.ratio.c$hmc_500_cut3 <- with(d0.comb.ratio.c, cut(hmc_500, c(min(hmc_500), 
    -1, 1, max(hmc_500))))
d0.comb.ratio$hmc_500_cut3 <- d0.comb.ratio.c$hmc_500_cut3[match(d0.comb.ratio$id, 
    d0.comb.ratio.c[, 1])]
```



```r
gg <- ggplot(na.omit(d0.comb.ratio[d0.comb.ratio$mod == "nuc" & d0.comb.ratio$pos == 
    "500", ]), aes(hmc_500_cut3, ko.wt))
gg + geom_boxplot() + coord_cartesian(ylim = c(-1, 2))
```

![plot of chunk unnamed-chunk-11](figure/unnamed-chunk-11.png) 



```r
with(d0.comb.ratio[d0.comb.ratio$mod == "nuc" & d0.comb.ratio$pos == "500", 
    ], wilcox.test(ko.wt[as.numeric(hmc_500_cut3) == 1], ko.wt[as.numeric(hmc_500_cut3) == 
    2]))
```

```
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  ko.wt[as.numeric(hmc_500_cut3) == 1] and ko.wt[as.numeric(hmc_500_cut3) == 2] 
## W = 212279, p-value = 0.00004158
## alternative hypothesis: true location shift is not equal to 0
```

