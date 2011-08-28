plot.path <- "~/storage/analysis/mprofiles/plots"
mp.path <- "~/storage/analysis/mprofiles"
anno.path <- "~/lib/annotations"
anno_hires.path <- "~/lib/annotations_hires"
feature.path <- "~/lib/features_general"
group2.path <- "~/storage/analysis/rna/quantiles"
profile2.path <- "~/s2/analysis/profiles"
medip_split.path <- "~/s2/data/medips_split"
group.classes <- c('character', rep('numeric', times=18), 'character', 'numeric')
gene.profiles <- c("gene_whole_W200N50F50", "transcSS_W200N50", "transcES_W200N50")
mp.classes <- c('character', 'character', rep('numeric', times=18))
profile.classes <- c('character', 'numeric', 'numeric', 'character', 'numeric', 'character', 'numeric', 'numeric')
samples.cells <- paste(c("omp", "ngn", "icam"), rep(c("hmedip.bed", "medip.bed"), each=3), sep="_")
samples <- c("moe_wt_hmc.bed", "moe_dnmt3a_hmc.bed")
samples.sc <- paste(c("cv", "iv", "cd", "id"), rep(c("hmc.bed", "mc.bed"), each=4), sep="_")
samples.d3a <- c("moe_wt_mc.bed", "moe_d3a_mc.bed", "moe_wt_hmc.bed", "moe_d3a_hmc.bed")
samples.cells.rlm <- paste(c("omp", "ngn", "icam"), rep(c("hmc_rlm", "mc_rlm"), each=3), sep="_")
samples.d3a.rlm <- paste(c("moe_wt", "moe_d3a"), rep(c("hmc_rlm", "mc_rlm"), each=2), sep="_")
