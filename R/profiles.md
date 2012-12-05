Plotting metaprofiles
========================================================

This is an R Markdown document. Markdown is a simple formatting syntax for authoring web pages (click the **MD** toolbar button for help on Markdown).

When you click the **Knit HTML** button a web page will be generated that includes both content as well as the output of any embedded R code chunks within the document. You can embed an R code chunk like this:


```r
suppressPackageStartupMessages(source("~/src/seqAnalysis/R/profiles2.R"))
```




```r
plot2("gene_whole_W200N50F50_chr", "omp_hmc_120424_rpkm", data_type = "rpkm/mean")
```

```
## [1] "omp_hmc_120424_rpkm"
## [1] "/media/storage2/analysis/profiles"
```

```
## [1] -0.193  0.171
```

![plot of chunk unnamed-chunk-2](figure/unnamed-chunk-2.png) 


