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
gg + geom_pointrange(aes(color=variable)) +
  facet_grid(ip~celltype) + opts(strip.background=theme_blank(),
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
