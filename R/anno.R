#####
# Annotation object containing:
#     BED positions
#     Names
#     Group
#     Strand
#     Coupling vector for a given nucleotide pattern
#
#####
library(utils)
library(MEDIPS)

source("~/src/R/AnnoSet-class.R")
source("~/src/R/util.R")

AnnoSet <- function(anno.name, set_chr=NULL, set_pos=NULL, set_names=NULL,
                    set_group=NULL, set_strand=NULL, chr_names=NULL, bin_size=NULL) {
  obj <- list(name = anno_name,
              set_chr = set_chr,
              set_pos = set_pos,
              set_names = set_names,
              set_group = set_group,
              set_strand = set_strand,
              chr_names = chr_names,
              bin_size = bin_size
              )
  class(obj) <- "AnnoSet"
  return(obj)
}

AnnoBedToAnnoData <- function(anno.bed.name) {
  anno.bed <- read.delim(anno.bed.name, header = FALSE,
                         col.names = c("chr", "start", "end", "name", "group", "strand"),
                         colClasses = c("character", "integer", "integer", "character",
                           "integer", "character"))
  chr.names <- sort(unique(anno.bed$chr))
  bin.size <- with(anno.bed, end[1] - start[1] + 1)
  pos <- anno.bed$start + 1

  AnnoSetObj <- new('AnnoSet', name = anno.bed.name, set_chr = anno.bed$chr, set_pos = pos,
                    set_names = anno.bed$name, set_group = anno.bed$group,
                    set_strand = anno.bed$strand, chr_names = chr.names, bin_size = bin.size)
  return(AnnoSetObj)
}

addDataIntersection <- function(annoSet, data.path) {
  data <- read.delim(data.path, header = FALSE, col.names = c("chr", "start", "end", "strand"),
                     colClasses = c("character", "integer", "integer", "character"))
  data.path.split <- strsplit(data.path,"/")
  data.path.name <- data.path.split[[1]][length(data.path.split[[1]])]
  counts <- MEDIPS.distributeReads(reads_start = data$start, reads_stop = data$end,
                                   reads_strand = data$strand, positions = annoSet@set_pos)
  annoSet <- addDataSet(annoSet, data.path.name, counts)
  return(annoSet)
}

generateProfile <- function(annoSet, sample=NULL, data.type="raw", FUN=mean) {
  seq_data <- annoSet@seq_data[[sample]]
  data_vec <- seq_data@data.type
  group <- annoSet@set_group
  result <- tapply(data_vec, group, FUN)
  return(result)
}

addData <- function(annoSet=NULL, data.name=NULL, counts=NULL) {
  annoSetgenome_data[[data.name]]['raw'] <- counts
  return(annoSet)
 
}

addPattern <- function(annoSet=NULL, pattern="CG", chr=NULL, start=NULL) {
   annoSet$pattern_data[[pattern]] <- cbind(chr = chr, start = start)
   return(annoSet)
 }
  
makeCouplingVector <- function(annoSet=NULL, pattern="CG", fragmentLength=700, func="count") {
  # pass AnnoSet param to MEDIPS.couplingVector(pattern, "count")
  setCoup <- AnnoSet.couplingVector(data = annoSet, pattern = pattern, fragmentLength = fragmentLength,
                         func = func )
  setCoup_name <- paste(fragmentLength, func, sep = "_")
  annoSet$coupling_vectors[[pattern]][setCoup_name] <- setCoup
  return(annoSet)
}

normalizeData <- function(AnnoSet=NULL, data, pattern="CG") {
  # normalize specified data by given coupling vector
}

  
