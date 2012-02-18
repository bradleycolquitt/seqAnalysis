library(Biostrings)
library(BSgenome)
library(BSgenome.Mmusculus.UCSC.mm9)

filterByLength <- function(seq, len) {
  return(seq[width(seq) >= len])
}

filterByNGylco <- function(seq) {
  pattern <- "N!S"  
  pattern_ind <- vmatchPattern(pattern, seq)
  return(pattern_ind)

}
