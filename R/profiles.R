####
# make profile of seq data signal over aligned features
# use ".profiles"
# tapply based on group
####

source("~/src/R/LEA/dipdata.R")
source("~/src/R/util.R")

#annotations_lib <- list.files("~/lib/annotations_old")
annotation_lib <- list.files("~/lib/annoset_export")
                                        #annotation_lib <- c("exon_ends_W200N25","exon_starts_W200N25","transcSS_W200N50","gene_whole_W200N50F50","transcES_W200N50")

hmedips <- paste(c("omp","ngn","icam"), "hmedip", sep="_")
medips <- paste(c("omp", "ngn", "icam"), "medip", sep="_")

labelDecode <- function(input) {
  if(input == "gene") return("gene_whole_W200N50F50")
  else if (input == "transcSS") return("transcSS_W200N50")
  else if (input == "transcES") return("transcES_W200N50")
  else if (input == "exon_starts") return("exon_starts_W200N50")
  else if (input == "exon_ends") return("exon_ends_W200N50")
}

intersectDipDataAnno <- function(dipdata.name, type, anno.name, FUN) {
  dipdata <- load.DipData(dipdata.name)
  genome <- dipdata$genome.data
  #anno.path <- "/home/user/lib/annotations"
  anno.path <- "/home/user/lib/annoset_export"
  anno <- read.delim(paste(anno.path, anno.name, sep="/"), header=FALSE,
                     colClasses=c('character', 'integer', 'integer', 'character', 'integer',
                                                'character'),
                     col.names=c('chr','start','end','name','group','strand')
                     )
  chrM.vec <- c("chrM", 0, 0, "", 0, "")
  anno <- rbind(anno, chrM.vec)

  genome.split <- bigsplit(genome, 'chr', splitcol = 2)
  anno.split <- split(anno, anno$chr)
  combined.group.index <- foreach(curr.genome= genome.split, curr.anno= anno.split,
                                  .combine = c) %dopar%  {
#    index <- rep(NA, times=length(curr.genome))
    index <- curr.anno[match(curr.genome, curr.anno$start),'group']
#    output.mat <- matrix(NA, dim=c(length(curr.genome),2))
#    output.vec <- rep(NA, times=length(curr.genome))
#    output.mat[,1] <- curr.anno[index, 'group']
    return(index)
  }
  #return(combined.group.index)
  result <- tapply(genome[,type], combined.group.index, mean)
  plotProfiles(result, paste(dipdata.name, type, anno.name, sep="_"))
  return(result)
#  result <- tapply(combined.groups[,2], combined.groups[,1], FUN)
}

generateProfiles <- function(samples, type, FUN=mean, write=FALSE) {
  output <- list()
  annotation_lib <- "tss_W200N100_export"
  for (sample in samples) {
    for (anno in annotation_lib) {
      #cat(paste("Working on ", paste(sample, type, sep="_"), 
      input.path <- paste("~/analysis/profiles", paste(sample, type, "profiles", sep="_"),
                          anno, sep = "/") 
      input <- read.delim(input.path, header = F,
                          col.names = c("chr", "start", "end",
                                        "name", "group", "strand", "reads"),
                          colClasses = c("character", "integer", "integer", "character",
                            "integer", "character", "integer"))
      result <- with(input, tapply(reads, group, FUN))
      if (write) {
        plotProfiles(result, paste(sample, type, paste("profiles/", anno, sep = ""), sep = "_"))
      } else {
        plotProfiles(result, title=anno)
      }  
      output <- lappend(output, c(result))
    }
  }  
  return(output)
}

plotProfiles <- function(profile, title=NULL, file_name=NULL, ...) {
  if (is.null(file_name)) {
    x11("", 7, 7)
  } else {
    file.path <- paste("~/analysis/profiles/", file_name, sep = "")
    png(paste(file.path, ".png", sep = ""))
  }  
  xlim = c(0,length(profile))
  ylim = c(floor(min(profile)), ceiling(max(profile)))
  plot(1, 1, type = "n", xlim = xlim, ylim = ylim, main=title )
  lines(profile, ...)
  if (!is.null(file_name)) dev.off()
  
}
