library(foreach)
library(doMC)

#plot.path <- "~/storage/analysis/mprofiles/plots"
plot.path <- "~/storage/analysis/mprofiles/features/plots"
mp.path <- "~/storage/analysis/mprofiles/features/norm"
anno.path <- "~/lib/annotations"
anno_hires.path <- "~/lib/annotations_hires"
feature.path <- "~/lib/features_general"
cluster.path <- "~/s2/analysis/cluster"
group2.path <- "~/storage/analysis/group2"
profile2.path <- "~/s2/analysis/profiles"
medip_split.path <- "~/s2/data/medips_split"
phase.path <- "~/s2/data/phasogram"
group.classes <- c('character', rep('numeric', times=18), 'character', 'numeric')
gene.profiles <- c("gene_whole_W200N50F50", "transcSS_W200N50", "transcES_W200N50")
mp.classes <- c('character', 'character', rep('numeric', times=18))
profile.classes <- c('character', 'numeric', 'numeric', 'character', 'numeric', 'character', 'numeric', 'numeric')
image.classes <- c('character', 'numeric', 'numeric', 'character', 'numeric', 'character', 'numeric')
samples.hmc <- c("omp_hmedip_ngn_hmedip", "omp_hmedip_icam_hmedip", "ngn_hmedip_icam_hmedip")
samples.mc <- c("omp_medip_ngn_medip", "omp_medip_icam_medip", "ngn_medip_icam_medip")
samples.cells <- paste(c("omp", "ngn", "icam"), rep(c("hmedip.bed", "medip.bed"), each=3), sep="_")
samples.cells_norm <- paste(c("omp", "ngn", "icam"), rep(c("hmc", "mc"), each=3), sep="_")
samples.cells_rpkm <- paste(c("omp", "ngn", "icam"), rep(c("hmc_rpkm", "mc_rpkm"), each=3), sep="_")
samples.cells_rpkm[1] <- "omp_hmc_120424_rpkm"
samples.cells_noM <- paste(c("omp", "ngn", "icam"), rep(c("hmc_noM", "mc_noM"), each=3), sep="_")
samples <- c("moe_wt_hmc.bed", "moe_dnmt3a_hmc.bed")
samples.ac3 <- c("moe_hmedip_moe_ac3_hmedip")
samples.sc <- paste(c("cv", "iv", "cd", "id"), rep(c("hmc.bed", "mc.bed"), each=4), sep="_")
samples.d3a <- c("moe_wt_mc.bed", "moe_d3a_mc.bed", "moe_wt_hmc.bed", "moe_d3a_hmc.bed")
samples.d3a_norm <- c("d3a_wt_hmc", "d3a_ko_hmc", "d3a_wt_mc", "d3a_ko_mc")
samples.d3a_rpkm <- c("moe_d3a_wt_hmc_rpkm", "moe_d3a_wt_hmc_30M_rpkm", "moe_d3a_ko_hmc_rpkm", "moe_d3a_wt_mc_rpkm", "moe_d3a_ko_mc_rpkm")
samples.d3a_2 <- c("moe_d3a_wt_hmc", "moe_d3a_ko_hmc", "moe_d3a_wt_mc", "moe_d3a_ko_mc")
samples.cells.rlm <- paste(c("omp", "ngn", "icam"), rep(c("hmc_rlm", "mc_rlm"), each=3), sep="_")
samples.d3a.rlm <- paste(c("moe_wt", "moe_d3a"), rep(c("hmc_rlm", "mc_rlm"), each=2), sep="_")
samples.d3a_nuc <- c("d3xog_wt_nuc_478_rmdup", "d3xog_ko_nuc_256_rmdup")
samples.d3a_nuc_sub <- c("d3xog_wt_nuc_478_rmdup_sub_omp_naked_s78", "d3xog_ko_nuc_256_rmdup_sub_omp_naked_s78")
samples.cells_raw <- c("omp_hmc_unnorm", "ngn_hmc_unnorm", "icam_hmc_unnorm",
                       "omp_mc_unnorm", "ngn_mc_unnorm", "icam_mc_unnorm")
samples.tfo_unnorm <- c("tfo_hmc_unnorm", "omp_hmc", "tfo_mc_unnorm", "omp_mc")
samples.tfo <- c("tfo_hmc", "omp_hmc", "tfo_mc", "omp_mc")
samples.medips_rf <- paste(c("omp", "ngn", "icam"), rep(c("hmc_rf", "mc_rf"), each=3), sep="_")
samples.nuc <- c("omp_nuc_0123", "icam_nuc_01234")
col4 <- brewer.pal(4, "Spectral")
col4_mod <- c(col4[1], brewer.pal(5, "PuOr")[1], brewer.pal(5, "RdYlGn")[5], col4[4])
col3 <- brewer.pal(3, "Dark2")
col2 <- rev(brewer.pal(3, "Set1")[1:2])

gene.features <- c("refgene_10to20kb_upstream", "refgene_2to10kb_upstream",
                   "refgene_1to3kb_up", "Refgene_1kb_up", "cgi", "Refgene_5_UTR",
                   "Refgene_exons",
                   "Refgene_intron", "Refgene_CDS", "Refgene_3_UTR", "Refgene_1kb_down",
                   "refgene_2to10kb_downstream", "refgene_10to20kb_downstream",
                   "intergenic_sub_rmsk")
moe.features <- c(gene.features, "MACS_DNaseHS_bed4.bed", "ngn_mk4_intergenic_inter_cons.bed4",
                  "omp_mk4_intergenic_inter_cons.bed4")
rmsk.features <- c("rmsk_DNA", "rmsk_LINE", "rmsk_LTR", "rmsk_rRNA",
                   "rmsk_SINE", "rmsk_snRNA")
cell.names <- c("icam", "ngn", "omp")
dip.names <- c("5hmC", "5mC")
