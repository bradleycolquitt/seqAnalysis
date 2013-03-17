O/Tet3 vs. Ngn
========================================================

* What genes are differentially expressed OMP vs. O/Tet3?
* What percentage of these correspond to OMP vs. Ngn DE?
* 





```r
suppressPackageStartupMessages(source("~/src/seqAnalysis/R/image.R"))
omp.hmc.pos <- makeImage("omp_hmc_120424_rpkm", "gene_whole_W200N50F50_chr", 
    data_type = "rpkm/mean", image = FALSE)
```

```
## [1] "/media/storage2/analysis/profiles/norm/rpkm/mean/gene_whole_W200N50F50_chr/images/omp_hmc_120424_rpkm"
```

```r
ott3.hmc.pos <- makeImage("ott3_1_hmc_rpkm", "gene_whole_W200N50F50_chr", data_type = "rpkm/mean", 
    image = FALSE)
```

```
## [1] "/media/storage2/analysis/profiles/norm/rpkm/mean/gene_whole_W200N50F50_chr/images/ott3_1_hmc_rpkm"
```

```r
ngn.hmc.pos <- makeImage("ngn_hmc_rpkm", "gene_whole_W200N50F50_chr", data_type = "rpkm/mean", 
    image = FALSE)
```

```
## [1] "/media/storage2/analysis/profiles/norm/rpkm/mean/gene_whole_W200N50F50_chr/images/ngn_hmc_rpkm"
```

```r
ngn.omp.hmc.pos <- ngn.hmc.pos - omp.hmc.pos
```




Principal components analysis on gene body to order

```r
omp.hmc.pos.pc <- prcomp(omp.hmc.pos[, 51:100])
ott3.hmc.pos.pc <- prcomp(ott3.hmc.pos[, 51:100])
omp.hmc.pos.pred <- predict(omp.hmc.pos.pc, omp.hmc.pos[, 51:100])
ott3.hmc.pos.pred <- predict(ott3.hmc.pos.pc, ott3.hmc.pos[, 51:100])

omp.hmc.pos.pc1 <- omp.hmc.pos[match(rownames(omp.hmc.pos.pred[order(omp.hmc.pos.pred[, 
    1]), ]), rownames(omp.hmc.pos)), ]
ott3.hmc.pos.pc1 <- ott3.hmc.pos[match(rownames(ott3.hmc.pos.pred[order(ott3.hmc.pos.pred[, 
    1]), ]), rownames(ott3.hmc.pos)), ]
ngn.hmc.pos.omp.pc1 <- ngn.hmc.pos[match(rownames(omp.hmc.pos.pred[order(omp.hmc.pos.pred[, 
    1]), ]), rownames(ngn.hmc.pos)), ]
omp.hmc.pos.ott3.pc1 <- omp.hmc.pos[match(rownames(ott3.hmc.pos.pred[order(ott3.hmc.pos.pred[, 
    1]), ]), rownames(omp.hmc.pos)), ]
ott3.hmc.pos.omp.pc1 <- ott3.hmc.pos[match(rownames(omp.hmc.pos.pred[order(omp.hmc.pos.pred[, 
    1]), ]), rownames(ott3.hmc.pos)), ]
```



```r
ngn.omp.hmc.pos.omp.pc1 <- na.omit(ngn.omp.hmc.pos[match(rownames(omp.hmc.pos.pc1), 
    rownames(ngn.omp.hmc.pos)), ])
```


Order difference matrix

```r
ott3.omp.hmc.pos.omp.pc1 <- na.omit(ott3.omp.hmc.pos[match(rownames(omp.hmc.pos.pc1), 
    rownames(ott3.omp.hmc.pos)), ])
```

```
## Error: object 'ott3.omp.hmc.pos' not found
```

```r
ott3.omp.hmc.pos.ott3.pc1 <- na.omit(ott3.omp.hmc.pos[match(rownames(ott3.hmc.pos.pc1), 
    rownames(ott3.omp.hmc.pos)), ])
```

```
## Error: object 'ott3.omp.hmc.pos' not found
```




```r
ngn.omp.hmc.pos.omp.pc1.c100 <- foreach(c = isplitRows(ngn.omp.hmc.pos.omp.pc1[, 
    51:100], chunks = 100), .combine = "rbind") %do% mean(c, na.rm = TRUE)
ngn.omp.hmc.pos.omp.pc1.c100 <- as.data.frame(ngn.omp.hmc.pos.omp.pc1.c100)
ngn.omp.hmc.pos.omp.pc1.c100$index <- 100:1
```



```r
MP.heat(ngn.hmc.pos.omp.pc1, range = c(0, 1.5), average = 50)
```

```
## Warning: executing %dopar% sequentially: no parallel backend registered
```

![plot of chunk gene_whole_W200N50F50_omp_hmc_120424_rpkm_mean_ordered_by_gene_body_pc1](figure/gene_whole_W200N50F50_omp_hmc_120424_rpkm_mean_ordered_by_gene_body_pc1.png) 


O/TT3 - OMP gene body 5hmC ordered by OMP 5hmC PC1

```r
gg <- ggplot(ngn.omp.hmc.pos.omp.pc1.c100, aes(V1, index))
```

```
## Error: could not find function "ggplot"
```

```r
gg + geom_point() + xlab("RPM") + ylab("") + theme(axis.text.y = element_blank()) + 
    coord_cartesian(xlim = c(-2, 0.5))
```

```
## Error: object 'gg' not found
```

