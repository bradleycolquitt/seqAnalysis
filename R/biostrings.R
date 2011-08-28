library(Biostrings)
library(BSgenome)
library(BSgenome.Mmusculus.UCSC.mm9)
library(foreach)
library(multicore)
library(doMC)

registerDoMC(cores=12)

chrs <- paste("chr", c(1:19, "X", "Y"), sep="")
BS.matchPattern <- function(positions, pattern) {
  positions.split <- split(positions, positions[,1])
  positions.split <- positions.split[-grep("random", names(positions.split))]
  seq <- foreach(chr=chrs, position=positions.split) %dopar% {
    curr.seq <- Views(Mmusculus[[chr]], start=position[,2], end=position[,3])
    dss <- DNAStringSet(curr.seq)
    names(dss) <- position[,4]
    return(dss)
  }
  
  counts <- foreach(curr.seq=seq) %dopar% {
    m <- vmatchPattern(pattern, curr.seq)
    counts <- countIndex(m)
    names(counts) <- names(curr.seq)
    return(counts)
  }
  return(counts)
}

