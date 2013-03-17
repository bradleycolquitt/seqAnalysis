library(ggplot2)

stat_sum_df <- function(fun, geom="crossbar", ...) {
  stat_summary(fun.data=fun, color="black", geom=geom, width=0.5, ...)
}

## For 5hmC, 5mC, mRNA associations

ggplot.2 <- function(data, label_thresh) {
  #ind_1 <- data$Value[data$Measure == "5hmC"] <= label_thresh & data$Value[data$Measure == "5hmC"] <= label_thresh
  #return(ind_1)
  #ind_2 <- data$Value[data$Measure == "5mC"] <= label_thresh
  size=3
  hjust=0
  vjust=-1
  gg <- ggplot(data, aes(Value[Measure=="5hmC"], Value[Measure=="5mC"]))
  gg +
  geom_point(alpha=I(1/3)) +
  opts(axis.title.x=theme_blank(), axis.title.y=theme_blank())
  #geom_text(aes(x=Value[Measure=="5hmC" & ind_1], y=Value[Measure=="mRNA" & ind_1], label=Gene[ind_1]), hjust=hjust, vjust=vjust, size=size, color="red")
  #+ geom_text(aes(x=Value[Measure=="5hmC" & ind_2], y=Value[Measure=="mRNA" & ind_2], label=Gene[ind_2]), hjust=.25, vjust=.25, size=size, color="blue") 
}

ggplot.box <- function(data, thresh) {
  data$thresh <- data$Value[data$Measure=="5hmC"] <= thresh
  gg <- ggplot(data, aes(thresh, Value[Measure=="mRNA"]))
  gg + geom_boxplot(aes(fill=thresh))
}

gg.scatter <- function(data, x, y, alpha, xlim=NULL, ylim=NULL) {
  ggplot(data, aes(x,y)) + geom_point(alpha=I(1/alpha)) 
}

#gg.filt + geom <- point(aes(y=Value[Measure=="5mC"]), alpha=I(1/5)) + geom <- text(aes(x=Value[Measure=="5hmC" & ind.hmc], y=Value[Measure=="5mC" & ind.hmc], label=Gene[ind.hmc], color="blue", size=.5, hjust=-.25, vjust=-.5))

#gg.rna + geom <- density(aes(x=FPKM, y=..density..)) + facet <- grid(var1000~.) + scale <- x <- continuous(limit=c(-10, 12))


#FOR RNA tracking plots
#rna.melt is ~/s2/analysis/rna/cells_mrna_mclust_4_and_5_melt
#gg.rna <- ggplot(rna.melt[!is.na(rna.melt$var2000.cl),], aes(variable, value, group=id, color=var2000.cl4))
#gg.rna + geom_line(alpha=I(1/10)) + facet_grid(var2000.cl4~.)

#FOR DNA boxplots split by RNA tracking classes
#dna.bx is ~/s2/analysis/features/norm/cells_hmc_mc_norm_mrna_var2000_cl4
#gg.bx <- ggplot(dna.bx, aes(cell, value, fill=class))
#gg.bx + geom_boxplot() + facet_grid(class~mod) + ylim(0,.3)

#For TFO/OMP hmc/mc comparison, pointrange
#tfo_iqr is output of statSummary.allIQR(set="tfo", value_type="unnorm/mean", transf="sqrt")
#gg_tfo_iqr <- ggplot(tfo_iqr, aes(feature_fac, median, ymin=low, ymax=up, fill=celltype))
gg_point <- function(gg) {
gg + geom_pointrange(aes(color=sample)) +
  facet_grid(ip~geno) + opts(strip.background=theme_blank(),
                                   strip.text.x=theme_blank(), strip.text.y=theme_blank(),
                                   axis.text.x=theme_blank()) +
                                     scale_color_manual(values=c("forestgreen", "darkgreen", "firebrick3", "darkred"))}
#current pdf is "~/s2/analysis/features/plots/tfo_omp_hmc_mc_pointrange.pdf"

#For cells hmc/mc split by d3a ko/wt loss of 5hmC BF greater 10
#dna.melt is "~/s2/analysis/features/norm/unnorm/mean/summaries/cells_hmc_mc_conditioned_by_d3a_refgene_noclust10_chr_hmc_bf_gt10"
#gg_dna <- ggplot(dna.melt, aes(celltype, value, fill=d3a_down))
#gg_dna + geom_boxplot(aes(y=sqrt(value)), outlier.shape=NA) + facet_grid(mod~d3a_down) +  opts(strip.background=theme_blank(), strip.text.x=theme_blank(), strip.text.y=theme_blank(), axis.text.x=theme_blank()) + ylim(0, .8)
#pdf is "~/s2/analysis/features/plots/omp_ngn_icam_hmc_mc_sqrt_d3a_mod_loss_noloss_0_08.pdf"

#For d3a hmc/mc refgene densities
#dna.d3a.melt is "~/s2/analysis/features/norm/unnorm/mean/summaries/d3a_refgene_chr_melt"
#gg.d3a <- ggplot(dna.d3a.melt, aes(x=sqrt(value), color=geno))
#gg.d3a + geom_density() + facet_grid(.~mod) + xlim(0, .6) + scale_color_manual(values=c("darkred", "darkblue")) + opts(panel.grid.major=theme_blank(), strip.background=theme_blank(), strip.text.x=theme_blank(), strip.text.y=theme_blank(), axis.text.x=theme_blank(), axis.title.x=theme_blank(), axis.title.y=theme_blank(), legend.position="none")
#pdf is "~/s2/analysis/features/plots/d3a_hmc_mc_sqrt_total_refgene_densities.pdf"

#For FFT, hmc, mc omp/icam plots
#nuc.melt is "~/s2/analysis/nuc/phasograms/refgene_omp_ngn_icam_hmc_mc_unnorm_omp_icam_fft_max_omp_rna_order_chunk100_melt", values ordered by omp FPKMs in omp_ngn_icam_mrna_nozero_log2, grouped into 100 bins and averaged within group
gg.fft <- function(data) {
gg.nuc <- ggplot(data[data$celltype!="ngn",], aes(index, value, color=measure))
gg.nuc + geom_point(alpha=I(1/2), color="grey") + geom_line(aes(group=1), stat="smooth") + facet_grid(measure~celltype, scales="free_y") + opts(panel.grid.major=theme_blank(), strip.background=theme_blank(), strip.text.x=theme_blank(), strip.text.y=theme_blank(), axis.text.x=theme_blank(), axis.title.x=theme_blank(), axis.title.y=theme_blank(), legend.position="none")
}
#plot is ~/s2/analysis/nuc/plots/omp_icam_fft_hmc_mc_indexed_by_increasing_omp_FPKM.pdf

gg.box <- function(data) {
  theme_set(theme_bw())
  gg.box <- ggplot(drop.levels(data[as.character(data$mk4.hmc)!="MID" & as.character(data$cell)!="ngn" & as.character(data$mod)!="mc", ]), aes(mk4.hmc, value, fill=mk4.hmc))
  gg.box + geom_boxplot(outlier.shape=NA) + facet_grid(mod~cell) + ylim(0, .5) + opts(strip.background=theme_blank(), strip.text.x=theme_blank())
}

gg.box2 <- function(data) {
  theme_set(theme_bw())
  gg.box <- ggplot(drop.levels(data[as.character(data$mk4.hmc)!="MID" & as.character(data$variable)!="ngn", ]), aes(factor(mk4.hmc), log(value,2), fill=factor(mk4.hmc)))
  gg.box + geom_boxplot(outlier.shape=NA) + facet_grid(.~variable)  + opts(strip.background=theme_blank(), strip.text.x=theme_blank())
}

#For OMP/ICAM RNA dotplot
#RNA is ~/s2/analysis/rna/summaries/omp_ngn_icam_mrna_dup_nozero_log2
#gg <- ggplot(rna, aes(icam, omp))
#colorMan <- scale_color_manual(values=c("black", "#D55E00", "#0072B2"))
#gg + geom_point(aes(color=factor(diff)), alpha=I(1/5)) + colorMan + opts(legend.position="none")

#For feature barplots
#tfo.2 object at ~/s2/analysis/features/norm/unnorm/mean/summaries/tfo_feature_toplot_norm_intergenic_sub_rmsk_sqrt_obj
gg.feature_bar <- function(gg) {
#gg.tfo.2 <- ggplot(droplevels(tfo.2[tfo.2$ip=="hmc",]), aes(celltype, value.median, ymax=up, ymin=low, fill=class))
theme_set(theme_bw())
  gg + geom_bar(position="dodge") + geom_linerange(position="dodge")  + facet_grid(ip~features.nice) + featureScaleFill4 + opts(legend.position="none") + ylab("Median square-root RPKM") + xlab("Cell type")

}

gg.feature_crossbar <- function(gg) {
#gg.tfo.2 <- ggplot(droplevels(tfo.2[tfo.2$ip=="hmc",]), aes(celltype, value.median, ymax=up, ymin=low, fill=class))
theme_set(theme_bw())
  gg + geom_crossbar(position="dodge") + facet_grid(ip~feature) #+ featureScaleFill
}
#For TT3 refgene plot split by omp_hmc/ott3_2_hmc Bayes Factor GT/LT 3
#tt3.melt object at ~/s2/analysis/features/norm/rpkm/mean/summaries/tt3_omp_ngn_icam_hmc_mc_refgene_chr_tt3_omp_2_hmc_bf.Rdata
#gg.tt3 <- ggplot(droplevels(tt3.melt[tt3.melt$sample!="ott3.1",]), aes(sample, value.sqrt))
#gg.tt3 + geom_boxplot(aes(fill=omp.2.hmc), outlier.shape=NA) + facet_grid(ip~omp.2.hmc) + coord_cartesian(ylim=c(0, 2.2)) + featureScaleFill3

#For TT3 exon histograms
#ex.in.m.trim.nd at ~/s2/analysis/features/norm/unnorm/mean/bf/tt3_omp_hmc_2_hmc_refgene_noclust10_sqrt_bf10_exon_counts_nodup.Rdata
#gg.ex <- ggplot(droplevels(ex.in.m.trim.nd[ex.in.m.trim.nd$bf10!="GT",]), aes(value, y=..density..))
#gg.ex + geom_histogram(aes(fill=bf10), color="black", binwidth=1) + facet_grid(bf10~variable) + coord_cartesian(xlim=c(0, 50)) + scale_fill_manual(values=col2)
#pdf(file="~/s2/analysis/features/plots/tt3_omp_hmc_2_hmc_refgene_noclust10_sqrt_bf10_exon_counts_nodup.pdf", 5, 5)

#For feature-peak intersects
#cell is "~/s2/data/homer/peaks/intersections/Rdata/cells_hmc_pairwise_peaks_feature_intersects_F3.Rdata"
#theme_set(theme_bw())
#pdf(file="~/s2/analysis/features/plots/feature_peaks_intersects_cells_hmc_pairwise_F3.pdf", 6, 8)
#cell.gg <- ggplot(cell, aes(peak_set.fac, internal_norm, fill=feature.factor))
#cell.gg + geom_bar() + scale_fill_hue(l=40, labels=features_merge_toplot_short)  + opts(legend.title=theme_blank())

#For feature-peak intersects as ranked fractions
#tt3.inter2 is ~/s2/data/homer/peaks/intersections/Rdata/ott3_2_hmc_gc_input_omp_hmc_gc_F3_feature_intersections2.Rdata
#tt3.inter2.gg <- ggplot(tt3.inter2, aes(factor(rank), y=internal_norm))
#tt3.inter2.gg + geom_bar() + scale_x_discrete(labels=feature_name)

#For feature scatter plot with median point
#cds.gg <- ggplot(cds, aes(omp_hmc_120424_rpkm, o.tt3.2_hmc_rpkm))
#cds.gg + geom_point(alpha=I(1/40)) + coord_cartesian(xlim=c(0,4), ylim=c(0,4)) + geom_abline(intercept=0, slope=1, color="blue") + geom_point(data=cds.med, aes(omp_hmc_120424_rpkm, o.tt3.2_hmc_rpkm), color="red", size=4)

#For density plots of moe d3a wt/ko hmc/mc/mrna/nuc read counts over exons
#comb.tlen.norm.m is ~/s2/analysis/features/norm/pileup/sum/summaries/d3a_hmc_mc_mrna_nuc_pileup_Refgene_exons_split2_chr_sum_normByExonLength_normByTotalReadCount_melt.Rdata
#comb.gg <- ggplot(comb.tlen.norm.m, aes(value, color=geno))
#comb.gg + geom_density() + facet_grid(type~hmc.gt10, scales="free_y") + coord_cartesian(xlim=c(0, .2)) + labs(x="Normalized read counts", y="Density", colour="Genotype")

#For boxplots of cells mRNA split by MOE Dnmt3a WT/KO genes with bayes factor >=10
#cells.m is ~/s2/analysis/rna/summaries/omp_ngn_icam_mrna_dup_biasCorrect_plus1_log2_moe_d3a_wt_ko_hmc_rpkm_refgene_nodUp_extend5kb_sqrt_bf_gt10.Rdata
#cells.gg <- ggplot(cells.m, aes(variable, value))
#cells.gg + geom_boxplot(aes(fill=d3a.bf.gt10), outlier.shape=NA) + facet_grid(.~d3a.bf.gt10) + coord_cartesian(ylim=c(0, 9)) + ylab("log2(FPKM + 1)") + xlab("Cell type") + scale_x_discrete(labels=c("HBC", "GBC", "mOSN")) + opts(legend.position="none")

#For boxplots of MOE Dnmt3a WT/KO mRNA split by MOE Dnmt3a WT/KO genes with bayes factor >=10
#d3a.melt is ~/s2/analysis/rna/summaries/moe_d3a_wt_ko_mrna_1log2_moe_d3a_wt_ko_hmc_rpkm_refgene_nodUp_extend5kb_sqrt_bf_gt10.Rdata
#d3a.gg <- ggplot(d3a.melt, aes(variable, value))
#d3a.gg + geom_boxplot(aes(fill=bf.gt10), outlier.shape=NA) + facet_grid(.~bf.gt10) + coord_cartesian(ylim=c(0, 10)) + ylab("log2(FPKM + 1)") + xlab("Genotype") + scale_x_discrete(labels=c("WT", "KO")) + opts(legend.position="none")

#For barplots of DIP qpcr data
#d3.na2 is olfactome:Documents/Research/methylation/qpcr/012811/dips.Rdata
#gg <- ggplot(d3.na[d3.na$ip!="in",], aes(cell, ct.mean))
#gg %+% d3.na2[d3.na2$ip!="in",] + geom_bar(aes(y=ct.mean, fill=cell)) + geom_errorbar(aes(ymin=ct.min, ymax=ct.max), width=.2) + facet_grid(ip~primer, scales="free_y") + scale_fill_manual(values=rev(col3)) + opts(strip.text.y=theme_text(face="bold", size=10), legend.position="none") + xlab("Cell type") + ylab("Cts normalized by input")
#pdf(file="olfactome:Documents/Research/methylation/qpcr/012811/hmc_mc_noClk4.pdf", 18, 9)

#For barplots of FPKM data for genes associated with DIP qPCR
#cell.melt is olfactome:Documents/Research/methylation/qpcr/012811/omp_ngn_icam_mrna_dup_biasCorrect_plus1_log2_dip_primers.Rdata
#rna.gg <- ggplot(cell.melt, aes(variable, value))
#rna.gg + geom_bar(aes(fill=variable)) + facet_grid(.~primer) + scale_fill_manual(values=col3) + xlab("Cell type") + ylab("log2(FPKM + 1)") + opts(legend.position="none")
#pdf(file="olfactome:Documents/Research/methylation/qpcr/012811/rna.pdf", 8, 4)

#For hmc levels versus transcriptional index
#0to500bp comb.tss is ~/s2/analysis/profiles/norm/unnorm/mean/refGene_noRandom_order_outsides2_tss_W25F400_chr/omp_ngn_icam_hmc_FPKM_index_mean_trim05_chunk100.Rdata
#9to11kb comb.mid is ~/s2/analysis/profiles/norm/unnorm/mean/refGene_noRandom_order_outsides2_tss_W25F400_chr/omp_ngn_icam_hmc_FPKM_index_9kbTo11kb_mean_trim05_chunk100.Rdata
#comb.mid.gg <- ggplot(comb.mid, aes(index, value))
#c <- comb.mid.gg + facet_grid(cell~., scales="free_y")
#c + stat_smooth(aes(colour=cell), method="lm", se=FALSE) + geom_point(color="grey") + xlab("Transcriptional index") + opts(legend.position="none")
#for 0to500bp, stat_smooth(aes(colour=cell), se=FALSE, fullrange=FALSE, span=.8)
#pdf(file="omp_ngn_icam_hmc_0to500bp_ordered_by_fpkm_chunk100.pdf", 4, 8)
#pdf(file="~/s2/analysis/profiles/plots/omp_ngn_icam_hmc_9kbto11kb_ordered_by_fpkm_chunk100.pdf", 4, 8)

#For nuc str vs. hmc/mc scatter/contour plots
#g + geom_point(alpha=.1)+ stat_density2d(geom="polygon", color="blue", aes(alpha=..level..)) +coord_cartesian(ylim=c(-3,4), xlim=c(-3, 3)) + scale_alpha_continuous(limits=c(0, .3), breaks=seq(0, .3, .025)) + theme_bw() + xlab("OMP 5mC") + ylab("OMP nucleosome stringency") + opts(legend.position="none")

#For BF densities
#gg.bf + geom <- density(aes(value)) + facet <- grid(variable~.) + coord <- cartesian(xlim=c(-400, 300), ylim=c(0, .006)) + xlab("KO vs. WT Bayes factor") + ylab("Density") + opts(strip.text.y=theme <- text(size=12, face="bold"), panel.margin=unit(1, "cm"))

#For boxplot hmc, mc enhancer gene fpkm
#gg <- ggplot(dna4.3mn, aes(celltype, value, fill=factor(l2)))
#gg + geom <- boxplot(outlier.shape=NA) + facet <- grid(sampe~l2, scales="free_y") + opts(legend.position="none")
