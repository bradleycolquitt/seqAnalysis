##########
# Create AnnoSet object from given annotation positions
##########
library(gtools)
library(bigmemory)
library(bigtabulate)
library(BSgenome)
library(BSgenome.Mmusculus.UCSC.mm9)

source("~/src/R/AnnoSet/AnnoSet-class.R")
source("~/src/R/util.R")

AnnoSet.initialize <- function(input_bed_name=NULL, bin_size, BSgenome=NULL) {

  # input_bed contains:
  #     chr
  #     start
  #     end
  #     name
  #     group
  #     strand
  # Note: add 1 to starts

  if(!file.exists(input_bed_name)) {
    stop("File doesn't exist.")
  }
  
  cat(paste("Loading annotation file ", input_bed_name, "...\n", sep=""))
  input_bed <- read.delim(input_bed_name, header=FALSE,
                          colClasses=c('character', 'numeric', 'numeric', 'character', 'numeric', 'character'))
  #cat("Load done.\n")	
  ## Sort chromosomes
  #cat("Sorting\n")
  #print(dim(input_bed))
  chr_names = NULL 
  if(length(unique(input_bed[,1]))>1){chr_names = as.character(mixedsort(unique(input_bed[,1])))}
  if(length(unique(input_bed[,1]))==1){chr_names = as.character(unique(input_bed[,1]))}
  
  ## Get chromosome lengths for all chromosomes within data set.
  cat(paste("Loading chromosome lengths for ",BSgenome, "...\n", sep=""))
  dataset=get(ls(paste("package:", BSgenome, sep="")))
  chr_lengths=as.numeric(sapply(chr_names, function(x){as.numeric(length(dataset[[x]]))}))
  set_chr <- input_bed[,1]
  #set_chr <- chrVecToNum(set_chr)
  set_pos <- input_bed[,2] + 1
  #set_names <- geneToNum(input_bed[,4])
  set_names <- input_bed[,4]
  set_group <-  input_bed[,5]
  #set_strand <- strandToNum(input_bed[,6])
  set_strand <- input_bed[,6]
  chr_names <- unique(input_bed[,1])
  
  #anno_matrix <- cbind(set_chr, set_pos, set_names, set_group, set_strand)
  
  AnnoSetObj <- new('AnnoSet', name = input_bed_name,
                    set_chr = set_chr,
                    set_pos = set_pos,
                    set_names = set_names,
                    set_group = set_group,
                    set_strand = set_strand,
                    #anno_data = as.big.matrix(anno_matrix),
                    chr_names = chr_names,
                    chr_lengths = chr_lengths,
                    bin_size = bin_size,
                    genome_name = BSgenome)
  return(AnnoSetObj)
}


save.AnnoSet <- function(anno_set) {
  data.path <- paste(anno_set$name, "annoset", sep=".")
  if (file.exists(data.path)) {
    unlink(data.path, recursive=TRUE)
  }
  dir.create(data.path)
  write.table(anno_set[c('name', 'chr.names', 'chr.lengths','bin.size')], file=paste(data.path, anno_set$name, sep="/"), sep=",")
  .writeDipData(anno_set$anno_data, data.path, "data")
}


.writeData <- function(data.col, data.path, colname) {
  cat(data.path)
  write.big.matrix(dip.data.col, paste(data.path, paste(colname, "txt", sep="."), sep="/"),
                   col.names=TRUE)
}

.readData <- function(colname, dtype) {
  col.path <- paste(colname, "txt", sep=".")
  col.bin.path <- paste(colname, "bin", sep=".")
  return(read.big.matrix(col.path, header=TRUE, 
                         type=dtype, backingfile=col.bin.path,
                         backingpath=".",
                         descriptorfile=paste(colname, "desc",sep=".")))
}

load.Data <- function(data.name) {
  .doLoad <- function() {
    obj <- AnnoSet(data.name)
    if (!file.exists("data.bin") ||
        !file.exists("data.desc")) {
      obj$data <- .readData("data", "double")    
    }
    else {
      obj$data <- attach.big.matrix(dget("data.desc"))
    }
    chr.data <- read.table(data.name, header=TRUE, sep=",")
    obj$chr.names <- chr.data$chr.names
    obj$chr.lengths <- chr.data$chr.lengths
    obj$bin.size <- chr.data$bin.size
    return(obj)
  }
  # kind of a kludge but w/e
  old.wd <- getwd()
  data.path <- paste(data.name, "annoset", sep=".")
  if (!file.exists(data.path)) {
    stop("Dipdata file does not exist")
  }
  setwd(data.path) 
  tryCatch(obj <- .doLoad(), finally=setwd(old.wd))
  return(obj)
}



modify.AnnoSet <- function() {}
