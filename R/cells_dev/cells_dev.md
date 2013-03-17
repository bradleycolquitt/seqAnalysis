Developmental 5hmC vs. RNA 
========================================================

```r
rna.cells <- read.delim("~/s2/analysis/rna/summaries/omp_ngn_icam_mrna_dup_biasCorrect_plus1_log2")
dna <- read.delim("~/s2/analysis/features/norm/rpkm/mean/summaries/cells_rpkm2_refgene_chr_sqrt")
m <- match(rownames(dna), rownames(rna.cells))
comb <- cbind(dna, rna.cells[m, ])
comb <- na.omit(comb)
omp.k2 <- kmeans(comb[, 7], centers = 2)
plot(density(comb[ngn.k2$cluster == 1, 8]))
```

```
## Error: object 'ngn.k2' not found
```

```r
cor(comb[omp.k2$cluster == 1, ])
```

```
##                     omp_hmc_120424_rpkm ngn_hmc_rpkm icam_hmc_rpkm
## omp_hmc_120424_rpkm              1.0000       0.8754        0.7757
## ngn_hmc_rpkm                     0.8754       1.0000        0.8304
## icam_hmc_rpkm                    0.7757       0.8304        1.0000
## omp_mc_rpkm                      0.6745       0.6789        0.6680
## ngn_mc_rpkm                      0.6749       0.6993        0.6671
## icam_mc_rpkm                     0.6546       0.6765        0.6459
## omp                              0.2146       0.1450        0.1167
## ngn                              0.3070       0.2554        0.1829
## icam                             0.2699       0.2679        0.3480
##                     omp_mc_rpkm ngn_mc_rpkm icam_mc_rpkm    omp    ngn
## omp_hmc_120424_rpkm      0.6745      0.6749       0.6546 0.2146 0.3070
## ngn_hmc_rpkm             0.6789      0.6993       0.6765 0.1450 0.2554
## icam_hmc_rpkm            0.6680      0.6671       0.6459 0.1167 0.1829
## omp_mc_rpkm              1.0000      0.8444       0.8449 0.1145 0.1587
## ngn_mc_rpkm              0.8444      1.0000       0.8347 0.1172 0.1493
## icam_mc_rpkm             0.8449      0.8347       1.0000 0.1159 0.1372
## omp                      0.1145      0.1172       0.1159 1.0000 0.6824
## ngn                      0.1587      0.1493       0.1372 0.6824 1.0000
## icam                     0.2116      0.1905       0.1696 0.2862 0.3321
##                       icam
## omp_hmc_120424_rpkm 0.2699
## ngn_hmc_rpkm        0.2679
## icam_hmc_rpkm       0.3480
## omp_mc_rpkm         0.2116
## ngn_mc_rpkm         0.1905
## icam_mc_rpkm        0.1696
## omp                 0.2862
## ngn                 0.3321
## icam                1.0000
```

```r
cor.test(comb[omp.k2$cluster == 1, 1], comb[omp.k2$cluster == 1, 7])$p.value
```

```
## [1] 0
```

```r

ngn.k2 <- kmeans(comb[, 8], centers = 2)
plot(density(comb[ngn.k2$cluster == 1, 8]))
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-11.png) 

```r
cor(comb[ngn.k2$cluster == 1, ])
```

```
##                     omp_hmc_120424_rpkm ngn_hmc_rpkm icam_hmc_rpkm
## omp_hmc_120424_rpkm             1.00000      0.88376     0.7424629
## ngn_hmc_rpkm                    0.88376      1.00000     0.8031784
## icam_hmc_rpkm                   0.74246      0.80318     1.0000000
## omp_mc_rpkm                     0.50294      0.54196     0.6047960
## ngn_mc_rpkm                     0.55368      0.61151     0.6262928
## icam_mc_rpkm                    0.61542      0.64951     0.6407537
## omp                             0.15368      0.04576     0.0003633
## ngn                             0.18312      0.16542     0.0727638
## icam                            0.07352      0.05315     0.2033337
##                     omp_mc_rpkm ngn_mc_rpkm icam_mc_rpkm        omp
## omp_hmc_120424_rpkm     0.50294     0.55368      0.61542  0.1536845
## ngn_hmc_rpkm            0.54196     0.61151      0.64951  0.0457576
## icam_hmc_rpkm           0.60480     0.62629      0.64075  0.0003633
## omp_mc_rpkm             1.00000     0.89961      0.89872 -0.2202144
## ngn_mc_rpkm             0.89961     1.00000      0.90729 -0.1103386
## icam_mc_rpkm            0.89872     0.90729      1.00000 -0.0911686
## omp                    -0.22021    -0.11034     -0.09117  1.0000000
## ngn                    -0.10989    -0.03619     -0.02687  0.7045750
## icam                   -0.01942    -0.00992     -0.03032  0.3965310
##                          ngn     icam
## omp_hmc_120424_rpkm  0.18312  0.07352
## ngn_hmc_rpkm         0.16542  0.05315
## icam_hmc_rpkm        0.07276  0.20333
## omp_mc_rpkm         -0.10989 -0.01942
## ngn_mc_rpkm         -0.03619 -0.00992
## icam_mc_rpkm        -0.02687 -0.03032
## omp                  0.70457  0.39653
## ngn                  1.00000  0.46865
## icam                 0.46865  1.00000
```

```r
cor.test(comb[omp.k2$cluster == 1, 2], comb[omp.k2$cluster == 1, 8])$p.value
```

```
## [1] 0
```

```r

icam.k2 <- kmeans(comb[, 9], centers = 2)
plot(density(comb[icam.k2$cluster == 1, 8]))
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-12.png) 

```r
cor(comb[icam.k2$cluster == 1, ])
```

```
##                     omp_hmc_120424_rpkm ngn_hmc_rpkm icam_hmc_rpkm
## omp_hmc_120424_rpkm             1.00000      0.88060       0.68874
## ngn_hmc_rpkm                    0.88060      1.00000       0.80850
## icam_hmc_rpkm                   0.68874      0.80850       1.00000
## omp_mc_rpkm                     0.48268      0.57823       0.64324
## ngn_mc_rpkm                     0.54960      0.64172       0.64861
## icam_mc_rpkm                    0.59915      0.66382       0.65289
## omp                             0.27197      0.06822      -0.12368
## ngn                             0.29194      0.13348      -0.09304
## icam                            0.04227      0.06764       0.14497
##                     omp_mc_rpkm ngn_mc_rpkm icam_mc_rpkm      omp      ngn
## omp_hmc_120424_rpkm     0.48268     0.54960      0.59915  0.27197  0.29194
## ngn_hmc_rpkm            0.57823     0.64172      0.66382  0.06822  0.13348
## icam_hmc_rpkm           0.64324     0.64861      0.65289 -0.12368 -0.09304
## omp_mc_rpkm             1.00000     0.91592      0.91236 -0.23773 -0.15528
## ngn_mc_rpkm             0.91592     1.00000      0.92593 -0.11585 -0.04708
## icam_mc_rpkm            0.91236     0.92593      1.00000 -0.10015 -0.03671
## omp                    -0.23773    -0.11585     -0.10015  1.00000  0.88086
## ngn                    -0.15528    -0.04708     -0.03671  0.88086  1.00000
## icam                   -0.04079    -0.03010     -0.05852  0.20543  0.23139
##                         icam
## omp_hmc_120424_rpkm  0.04227
## ngn_hmc_rpkm         0.06764
## icam_hmc_rpkm        0.14497
## omp_mc_rpkm         -0.04079
## ngn_mc_rpkm         -0.03010
## icam_mc_rpkm        -0.05852
## omp                  0.20543
## ngn                  0.23139
## icam                 1.00000
```

```r
cor.test(comb[omp.k2$cluster == 1, 3], comb[omp.k2$cluster == 1, 9])$p.value
```

```
## [1] 0
```

```r

cor.test(comb.nzo[, 1], comb.nzo[, 7])$p.value
```

```
## Error: object 'comb.nzo' not found
```

```r
cor.test(comb.nzn[, 2], comb.nzn[, 8])$p.value
```

```
## Error: object 'comb.nzn' not found
```

```r
cor.test(comb.nzi[, 3], comb.nzi[, 9])$p.value
```

```
## Error: object 'comb.nzi' not found
```

