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


FOR RNA tracking plots
rna.melt is ~/s2/analysis/rna/cells_mrna_mclust_4_and_5_melt
gg.rna <- ggplot(rna.melt[!is.na(rna.melt$var2000.cl),], aes(variable, value, group=id, color=var2000.cl4))
gg.rna + geom_line(alpha=I(1/10)) + facet_grid(var2000.cl4~.)

FOR DNA boxplots split by RNA tracking classes
dna.bx is ~/s2/analysis/features/norm/cells_hmc_mc_norm_mrna_var2000_cl4
gg.bx <- ggplot(dna.bx, aes(cell, value, fill=class))
gg.bx + geom_boxplot() + facet_grid(class~mod) + ylim(0,.3)
