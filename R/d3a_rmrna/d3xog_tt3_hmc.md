D3xOG , O/Tet3 5hmC comparison
========================================================

Compare 5hmC levels in Dnmt3a WT/Het/KO and O/Tet3 over
* DNase HS peaks
* 2-7 kb downstream
* 1-3 kb upstream


```r
opts_chunk$set(warning = FALSE, message = FALSE, error = FALSE)
source("~/src/seqAnalysis/R/features.R")
```

```
## Loading required package: iterators
```

```
## Loading required package: multicore
```

```
## Loading required package: BiocGenerics
```

```
## Attaching package: 'BiocGenerics'
```

```
## The following object(s) are masked from 'package:stats':
## 
## xtabs
```

```
## The following object(s) are masked from 'package:base':
## 
## Filter, Find, Map, Position, Reduce, anyDuplicated, cbind, colnames,
## duplicated, eval, get, intersect, lapply, mapply, mget, order, paste,
## pmax, pmax.int, pmin, pmin.int, rbind, rep.int, rownames, sapply, setdiff,
## table, tapply, union, unique
```

```
## Loading required package: BSgenome
```

```
## Loading required package: IRanges
```

```
## Loading required package: GenomicRanges
```

```
## Loading required package: Biostrings
```

```
## Loading required package: Biobase
```

```
## Welcome to Bioconductor
## 
## Vignettes contain introductory material; view with 'browseVignettes()'. To
## cite Bioconductor, see 'citation("Biobase")', and for packages
## 'citation("pkgname")'.
```

```
## Loading required package: gtools
```

```
## Attaching package: 'gtools'
```

```
## The following object(s) are masked from 'package:boot':
## 
## inv.logit, logit
```

```
## Loading required package: gdata
```

```
## gdata: read.xls support for 'XLS' (Excel 97-2004) files ENABLED.
```

```
## ```

```
## gdata: read.xls support for 'XLSX' (Excel 2007+) files ENABLED.
```

```
## Attaching package: 'gdata'
```

```
## The following object(s) are masked from 'package:Biobase':
## 
## combine
```

```
## The following object(s) are masked from 'package:IRanges':
## 
## trim
```

```
## The following object(s) are masked from 'package:BiocGenerics':
## 
## combine
```

```
## The following object(s) are masked from 'package:stats':
## 
## nobs
```

```
## The following object(s) are masked from 'package:utils':
## 
## object.size
```

```
## Loading required package: caTools
```

```
## Attaching package: 'caTools'
```

```
## The following object(s) are masked from 'package:IRanges':
## 
## runmean
```

```
## Loading required package: grid
```

```
## Loading required package: KernSmooth
```

```
## KernSmooth 2.23 loaded Copyright M. P. Wand 1997-2009
```

```
## Loading required package: MASS
```

```
## Attaching package: 'gplots'
```

```
## The following object(s) are masked from 'package:IRanges':
## 
## space
```

```
## The following object(s) are masked from 'package:stats':
## 
## lowess
```

```
## Attaching package: 'plyr'
```

```
## The following object(s) are masked from 'package:IRanges':
## 
## compact, desc, rename
```

```
## Attaching package: 'reshape'
```

```
## The following object(s) are masked from 'package:plyr':
## 
## rename, round_any
```

```
## The following object(s) are masked from 'package:IRanges':
## 
## expand, rename
```

```
## Attaching package: 'proxy'
```

```
## The following object(s) are masked from 'package:stats':
## 
## as.dist, dist
```

```r
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
library(ggplot2)
```


DHS
--------

```r
hmc <- makeFeatureMatrix2("d3a_het_dnase_sort_q30_dnase_sub_igenome_ensembl_genes_extend5kb.bed_chr", 
    "d3xog_tt3_hmc", data_type = "rpkm/mean", select = "flank500")
```

```
##              omp_hmc_rep1_mean_omp_hmc_rep2 
##   "omp_hmc_rep1_mean_omp_hmc_rep2_flank500" 
##                    d3xog_het_hmc_paired_q30 
##         "d3xog_het_hmc_paired_q30_flank500" 
##                     d3xog_ko_hmc_paired_q30 
##          "d3xog_ko_hmc_paired_q30_flank500" 
##            ott3_hmc_rep1_mean_ott3_hmc_rep2 
## "ott3_hmc_rep1_mean_ott3_hmc_rep2_flank500"
```

```r
colnames(hmc) <- c("wt", "het", "ko", "tt3")
hmc.m <- melt(hmc)
hmc.m$id <- rownames(hmc)
```



```r
hmc.m <- ddply(hmc.m, .(id), mutate, value.norm.wt = (value + 0.01)/(value[variable == 
    "wt"] + 0.01))
```



```r
gg <- ggplot(hmc.m, aes(value, color = variable))
gg + geom_density() + scale_color_brewer(palette = "Set1") + coord_cartesian(xlim = c(0, 
    4))
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4.png) 



```r
gg <- ggplot(hmc.m, aes(variable, value.norm.wt, group = id))
gg + geom_line(alpha = I(1/10))
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5.png) 





```r
gg <- ggplot(hmc.m, aes(value.norm.wt, color = variable))
gg + geom_density() + facet_grid(variable ~ .) + coord_cartesian(xlim = c(0, 
    5))
```

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6.png) 



```r
gg <- ggplot(hmc.m, aes(variable, value.norm.wt))
gg + geom_boxplot() + coord_cartesian(ylim = c(0, 2))
```

![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-7.png) 


2 to 7 kb downstream 
-----------------------

```r
hmc.up <- makeFeatureMatrix2("refgene_nodup_TSS2to7kb.bed_chr", "d3xog_tt3_hmc", 
    data_type = "rpkm/mean")
```

```
## [1] "omp_hmc_rep1_mean_omp_hmc_rep2"   "d3xog_het_hmc_paired_q30"        
## [3] "d3xog_ko_hmc_paired_q30"          "ott3_hmc_rep1_mean_ott3_hmc_rep2"
```

```r
colnames(hmc.up) <- c("wt", "het", "ko", "tt3")
hmc.up.m <- melt(hmc.up)
hmc.up.m$id <- rownames(hmc.up)
```



```r
gg <- ggplot(hmc.up.m, aes(value, color = variable))
gg + geom_density() + scale_color_brewer(palette = "Set1") + coord_cartesian(xlim = c(0, 
    4))
```

![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-9.png) 


1 to 3 kb upstream
-----------------------

```r
hmc.gene <- makeFeatureMatrix2("refgene_1to3kb_up_chr", "d3xog_tt3_hmc", data_type = "rpkm/mean")
```

```
## [1] "omp_hmc_rep1_mean_omp_hmc_rep2"   "d3xog_het_hmc_paired_q30"        
## [3] "d3xog_ko_hmc_paired_q30"          "ott3_hmc_rep1_mean_ott3_hmc_rep2"
```

```r
colnames(hmc.gene) <- c("wt", "het", "ko", "tt3")
hmc.gene.m <- melt(hmc.gene)
hmc.gene.m$id <- rownames(hmc.gene)
```



```r
gg <- ggplot(hmc.gene.m, aes(value, color = variable))
gg + geom_density() + scale_color_brewer(palette = "Set1") + coord_cartesian(xlim = c(0, 
    4))
```

![plot of chunk unnamed-chunk-11](figure/unnamed-chunk-11.png) 

