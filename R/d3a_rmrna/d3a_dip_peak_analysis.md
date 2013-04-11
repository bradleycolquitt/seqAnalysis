D3a 5hmC 5mC peak analysis
========================================================

Peak intersections
-------------------


```r
suppressPackageStartupMessages(source("~/src/seqAnalysis/R/features.R"))
opts_knit$set(progress = TRUE, verbose = TRUE)
```



```r

hmc.wt.ko.feat <- processIntersectSummary("~/s2/data/homer/peaks/intersections/d3a/moe_d3a_wt_hmc_gc_input_moe_d3a_wt_in_gc_size1kb_F2_inter_moe_d3a_wt_hmc_gc_input_moe_d3a_ko_hmc_gc_size1kb_F2.bed/summary")
mc.wt.ko.feat <- processIntersectSummary("~/s2/data/homer/peaks/intersections/d3a/moe_d3a_wt_mc_gc_input_moe_d3a_wt_in_gc_size1kb_F2_inter_moe_d3a_wt_mc_gc_input_moe_d3a_ko_mc_gc_size1kb_F2.bed/summary")
# hmc.wt.ko.feat <- omp.ott3.feat[-grep('mOSN enhancer',
# omp.ott3.feat$feature.pretty),]
hmc.wt.ko.feat
```

```
##                                                                                     feature
## 1                                       phastCons30way_intergenic_merge500_thresh500_merged
## 2  phastCons30way_intergenic_sorted_merge500_thresh500_inter_omp_h3k4me1_default.bed_merged
## 3                                                                        Refgene_CDS_merged
## 4                                                                      Refgene_5_UTR_merged
## 5                                                                      Refgene_3_UTR_merged
## 6                                                                     Refgene_intron_merged
## 7                                                                  refgene_1to3kb_up_merged
## 8                                                                intergenic_sub_rmsk_merged
## 9                                                                     Refgene_1kb_up_merged
## 10                                                                               cgi_merged
##    reference_num feature_num feature_span count fraction fraction_norm
## 1           7041       58869     64778958   227 0.032240     0.0004977
## 2           7041        4094      4414009   106 0.015055     0.0034107
## 3           7041      189961     33521649    78 0.011078     0.0003305
## 4           7041       31208      5041097     9 0.001278     0.0002536
## 5           7041       21921     22959114    84 0.011930     0.0005196
## 6           7041      177062    881598741  4276 0.607300     0.0006889
## 7           7041       21900     44931696   364 0.051697     0.0011506
## 8           7041     2437128    905159578  1634 0.232069     0.0002564
## 9           7041       21334     22385987   143 0.020310     0.0009072
## 10          7041       16026     10496250    11 0.001562     0.0001488
##    internal_norm expected log2.obs.exp                          fisher.p
## 1        0.06096      229     -0.01266 0.5555648813933142804444287321530
## 2        0.41777       16      2.72792 0.0000000000000000121357980004362
## 3        0.04048      118     -0.59724 0.9982443247124013208093629145878
## 4        0.03106       18     -1.00000 0.9737912542561228956827790170792
## 5        0.06365       81      0.05247 0.4388803171745087605692958732106
## 6        0.08438     3111      0.45888 0.0000000000000000000000000001978
## 7        0.14093      159      1.19491 0.0000000000000000003404691717521
## 8        0.03140     3194     -0.96696 1.0000000000000000000000000000000
## 9        0.11113       79      0.85609 0.0000121022037138068422503075702
## 10       0.01823       37     -1.75002 0.9999683606600545671128088542901
##                          fisher.fdr
## 1  0.925941468988857208088916195265
## 2  0.000000000000000040452660001454
## 3  1.000000000000000000000000000000
## 4  1.000000000000000000000000000000
## 5  0.877760634349017521138591746421
## 6  0.000000000000000000000000001978
## 7  0.000000000000000001702345858761
## 8  1.000000000000000000000000000000
## 9  0.000030255509284517106472801873
## 10 1.000000000000000000000000000000
##                                                                              feature.factor
## 1                                       phastCons30way_intergenic_merge500_thresh500_merged
## 2  phastCons30way_intergenic_sorted_merge500_thresh500_inter_omp_h3k4me1_default.bed_merged
## 3                                                                        Refgene_CDS_merged
## 4                                                                      Refgene_5_UTR_merged
## 5                                                                      Refgene_3_UTR_merged
## 6                                                                     Refgene_intron_merged
## 7                                                                  refgene_1to3kb_up_merged
## 8                                                                intergenic_sub_rmsk_merged
## 9                                                                     Refgene_1kb_up_merged
## 10                                                                               cgi_merged
##          feature.pretty      class
## 1  Intergenic conserved Intergenic
## 2         mOSN enhancer Intergenic
## 3                   CDS      Exons
## 4                5' UTR      Exons
## 5                3' UTR      Exons
## 6                Intron Transcript
## 7    1 to 3 kb upstream   Upstream
## 8            Intergenic Intergenic
## 9          1kb upstream   Upstream
## 10                  CGI   Upstream
```

```r

mc.wt.ko.feat
```

```
##                                                                                     feature
## 1                                       phastCons30way_intergenic_merge500_thresh500_merged
## 2  phastCons30way_intergenic_sorted_merge500_thresh500_inter_omp_h3k4me1_default.bed_merged
## 3                                                                        Refgene_CDS_merged
## 4                                                                      Refgene_5_UTR_merged
## 5                                                                      Refgene_3_UTR_merged
## 6                                                                     Refgene_intron_merged
## 7                                                                  refgene_1to3kb_up_merged
## 8                                                                intergenic_sub_rmsk_merged
## 9                                                                     Refgene_1kb_up_merged
## 10                                                                               cgi_merged
##    reference_num feature_num feature_span count  fraction fraction_norm
## 1          10512       58869     64778958   233 0.0221651     0.0003422
## 2          10512        4094      4414009    21 0.0019977     0.0004526
## 3          10512      189961     33521649   410 0.0390030     0.0011635
## 4          10512       31208      5041097    10 0.0009513     0.0001887
## 5          10512       21921     22959114   203 0.0193113     0.0008411
## 6          10512      177062    881598741  4572 0.4349315     0.0004933
## 7          10512       21900     44931696   200 0.0190259     0.0004234
## 8          10512     2437128    905159578  1988 0.1891172     0.0002089
## 9          10512       21334     22385987    76 0.0072298     0.0003230
## 10         10512       16026     10496250    25 0.0023782     0.0002266
##    internal_norm expected log2.obs.exp                    fisher.p
## 1        0.07337      341     -0.54944 0.9999966527385157899487922
## 2        0.09705       23     -0.13124 0.6741491002454306213920177
## 3        0.24950      177      1.21187 0.0000000000000000000005035
## 4        0.04047       27     -1.43296 0.9987087326296759659527424
## 5        0.18037      121      0.74647 0.0000035975721214440935090
## 6        0.10579     4645     -0.02285 0.7410182793379864740757057
## 7        0.09080      237     -0.24489 0.9641461979078356225869584
## 8        0.04480     4769     -1.26237 1.0000000000000000000000000
## 9        0.06926      118     -0.63472 0.9989832255887917122905151
## 10       0.04859       55     -1.13750 0.9997690180349407818027885
##                    fisher.fdr
## 1  1.000000000000000000000000
## 2  1.000000000000000000000000
## 3  0.000000000000000000005035
## 4  1.000000000000000000000000
## 5  0.000017987860607220465851
## 6  1.000000000000000000000000
## 7  1.000000000000000000000000
## 8  1.000000000000000000000000
## 9  1.000000000000000000000000
## 10 1.000000000000000000000000
##                                                                              feature.factor
## 1                                       phastCons30way_intergenic_merge500_thresh500_merged
## 2  phastCons30way_intergenic_sorted_merge500_thresh500_inter_omp_h3k4me1_default.bed_merged
## 3                                                                        Refgene_CDS_merged
## 4                                                                      Refgene_5_UTR_merged
## 5                                                                      Refgene_3_UTR_merged
## 6                                                                     Refgene_intron_merged
## 7                                                                  refgene_1to3kb_up_merged
## 8                                                                intergenic_sub_rmsk_merged
## 9                                                                     Refgene_1kb_up_merged
## 10                                                                               cgi_merged
##          feature.pretty      class
## 1  Intergenic conserved Intergenic
## 2         mOSN enhancer Intergenic
## 3                   CDS      Exons
## 4                5' UTR      Exons
## 5                3' UTR      Exons
## 6                Intron Transcript
## 7    1 to 3 kb upstream   Upstream
## 8            Intergenic Intergenic
## 9          1kb upstream   Upstream
## 10                  CGI   Upstream
```



```r
hmc.wt.ko.feat$feature.pretty <- factor(hmc.wt.ko.feat$feature.pretty, levels = levels(hmc.wt.ko.feat$feature.pretty)[length(levels(hmc.wt.ko.feat$feature.pretty)):1])
mc.wt.ko.feat$feature.pretty <- factor(mc.wt.ko.feat$feature.pretty, levels = levels(mc.wt.ko.feat$feature.pretty)[length(levels(mc.wt.ko.feat$feature.pretty)):1])
```



```r
theme_set(theme_bw())

gg <- ggplot(hmc.wt.ko.feat, aes(feature.pretty, log2.obs.exp, fill = class))
gg + geom_bar(width = 0.8) + scale_fill_brewer(palette = "Set1") + theme(legend.position = "none") + 
    xlab("") + coord_flip() + ylab("log2 observed / expected") + geom_hline(yintercept = 0)
```

```
## Mapping a variable to y and also using stat="bin".  With stat="bin", it
## will attempt to set the y value to the count of cases in each group.  This
## can result in unexpected behavior and will not be allowed in a future
## version of ggplot2.  If you want y to represent counts of cases, use
## stat="bin" and don't map a variable to y.  If you want y to represent
## values in the data, use stat="identity".  See ?geom_bar for examples.
## (Deprecated; last used in version 0.9.2)
```

```
## Warning: Stacking not well defined when ymin != 0
```

![plot of chunk d3a_moe_wt_ko_hmc_feature_intersections_noMk4_obsExp_bar](figure/d3a_moe_wt_ko_hmc_feature_intersections_noMk4_obsExp_bar.png) 



```r
theme_set(theme_bw())

gg <- ggplot(mc.wt.ko.feat, aes(feature.pretty, log2.obs.exp, fill = class))
gg + geom_bar(width = 0.8) + scale_fill_brewer(palette = "Set1") + theme(legend.position = "none") + 
    xlab("") + coord_flip() + ylab("log2 observed / expected") + geom_hline(yintercept = 0)
```

```
## Mapping a variable to y and also using stat="bin".  With stat="bin", it
## will attempt to set the y value to the count of cases in each group.  This
## can result in unexpected behavior and will not be allowed in a future
## version of ggplot2.  If you want y to represent counts of cases, use
## stat="bin" and don't map a variable to y.  If you want y to represent
## values in the data, use stat="identity".  See ?geom_bar for examples.
## (Deprecated; last used in version 0.9.2)
```

```
## Warning: Stacking not well defined when ymin != 0
```

![plot of chunk d3a_moe_wt_ko_mc_feature_intersections_noMk4_obsExp_bar](figure/d3a_moe_wt_ko_mc_feature_intersections_noMk4_obsExp_bar.png) 

