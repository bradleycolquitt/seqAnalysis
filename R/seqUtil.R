library(BSgenome)
library(BSgenome.Mmusculus.UCSC.mm9)
library(foreach)

shiftBedPositions <- function(bed, shift, direction="up") {
  if (direction == "up") {
    direct <- 1
  } else {
    direct <- -1
  }
  left <- vector(length=nrow(bed))
  right <- vector(length=nrow(bed))
  plus.ind <- na.omit(bed[,6] == "+")
  #minus.ind <- bed[,6] == "-"
  #return(plus.ind)
  #print(shift)
  #return(bed[plus.ind, 2])
  left[plus.ind] <- bed[plus.ind, 2] - shift
  #return(left)
  right[plus.ind] <- bed[plus.ind, 2]
  right[!plus.ind] <- bed[!plus.ind, 3] + shift
  left[!plus.ind] <- bed[!plus.ind, 3]
  out <- cbind(bed[,1], left, right, bed[,4:6])
  return(out)
}
#take matrix with sequence names and sequence, and generate
#fasta file
formatFasta <- function(data, fname=NULL) {
  names <- data[,1]
  seq <- data[,2]
  fc <- file(fname, 'w')
  for (i in 1:nrow(data)) {
    name.out <- paste(c(">", names[i], "\n"), sep="")
    cat(name.out, file=fc)
    seq.out <- paste(c(seq[i], "\n"), sep="")
    cat(seq.out, file=fc)
  }
  close(fc)
}

trimBed <- function(bed, amt, pos="up") {
  plus.ind <- bed[,6] == "+"
  if (pos == "up") {
    plus.pos <- 2
    minus.pos <- 3
  } else {
    plus.pos <- 3
    minus.pos <- 2
  }
  bed[plus.ind, plus.pos] <- bed[plus.ind, plus.pos] + amt
  bed[!plus.ind, minus.pos] <- bed[!plus.ind, minus.pos] - amt
  return(bed)
}

countBasesByStrand <- function(seq, strand, pattern) {
  plus_ind <- strand == "+"
  up <- DNAString(pattern)
  down <- reverseComplement(up)
  counts <- list()
  for (i in 1:length(seq)) {
    if (plus_ind[i]) {
      up_count <- countPattern(up, seq[[i]])
      down_count <- countPattern(down, seq[[i]])
    } else {
      up_count <- countPattern(down, seq[[i]])
      down_cout <- countPattern(up, seq[[i]])
    }
    counts <- c(counts, list(c(up_count, down_count)))
  }
  return(counts)

}

freqBasesByStrand <- function(seq, strand, phase=0) {
  print(paste("Phase: ", phase, sep=""))
   plus_ind <- strand == "+"
  #up <- DNAString(pattern)
  #down <- reverseComplement(up)
  freqs <- list()
  pb <- txtProgressBar(min = 0, max = length(seq), style=3)
  for (i in 1:length(seq)) {
    setTxtProgressBar(pb, i)
    #print(nchar(seq[i]))
    if (nchar(seq[i]) > 10 & nchar(seq[i]) < 1200) {
      if (phase > 0) {
        phased_seq <- str_sub(seq[[i]], phase, phase)
        for (j in seq(phase + 3, nchar(seq[[i]]), by=3)) {
          phased_seq <- c(phased_seq, str_sub(seq[[i]], j, j))
        }
        phased_seq <- paste(phased_seq, collapse="")
        #print(phased_seq)
        seq_string <- DNAString(unlist(phased_seq))
      } else {
        seq_string <- DNAString(seq[[i]])
      }
      if (plus_ind[i]) {
        freq <- alphabetFrequency(seq_string, as.prob=TRUE)
      } else {
        freq <- alphabetFrequency(reverseComplement(seq_string), as.prob=TRUE)
      }
      freqs <- c(freqs, list(freq))
    }
  }
  return (ldply(freqs))
}

freqBasesByStrand.all3 <- function(seq, strand) {
  result <- foreach(phase=c(1:3)) %dopar% {
    
    return(freqBasesByStrand(seq, strand, phase))
  }
  return(result)
}



