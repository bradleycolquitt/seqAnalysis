library(BSgenome)
library(BSgenome.Mmusculus.UCSC.mm9)
library(foreach)
library(itertools)
library(doMC)

registerDoMC(cores=6)

bedClasses <- c("character", "numeric", "numeric", "character", "numeric", "character")
nuc_seq_path <- "/media/storage2/analysis/nuc/sequence"

shiftBedPositions <- function(bed, shift, pos="start", direction="up") {
  
  left <- vector(length=nrow(bed))
  right <- vector(length=nrow(bed))
  plus.ind <- bed[,6] == "+"

  if (pos=="start") {
    if (direction=="up") {
      left[plus.ind] <- bed[plus.ind, 2] - shift
      right[plus.ind] <- bed[plus.ind, 2]
      right[!plus.ind] <- bed[!plus.ind, 3] + shift
      left[!plus.ind] <- bed[!plus.ind, 3]
    } else if (direction=="down") {
      left[plus.ind] <- bed[plus.ind, 2]
      right[plus.ind] <- bed[plus.ind, 2] + shift
      right[!plus.ind] <- bed[!plus.ind, 3] 
      left[!plus.ind] <- bed[!plus.ind, 3] - shift
    }  
  } else if (pos=="end") {
    if (direction=="up") {
      left[plus.ind] <- bed[plus.ind, 3] - shift 
      right[plus.ind] <- bed[plus.ind, 3] 
      right[!plus.ind] <- bed[!plus.ind, 2] + shift
      left[!plus.ind] <- bed[!plus.ind, 2] 
    } else if (direction=="down") {
      left[plus.ind] <- bed[plus.ind, 3]
      right[plus.ind] <- bed[plus.ind, 3] + shift
      right[!plus.ind] <- bed[!plus.ind, 2] 
      left[!plus.ind] <- bed[!plus.ind, 2] - shift
    }  
  }
   
  out <- cbind(bed[,1], left, right, bed[,4:6])
  return(out)
}

#To adjust coordinates of all beds within a directory
shiftBedPositions.dir <- function(dir, shift, pos="start", direction="up") {
  files <- list.files(dir)
  ind <- grep("(start|end)", files)
  if (!is.null(ind)) files <- files[-ind]
  foreach(file=files) %dopar% {
    if (file=="subset") return
    bed <- read.delim(paste(dir, file, sep="/"), header=FALSE)
    bed <- shiftBedPositions(bed, shift=shift, pos=pos, direction=direction)
    write.table(bed, file=paste(dir, paste(file, pos, direction, shift, sep="_"), sep="/"),
                quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
  }
}

batchShift <- function(path, shift, pos="start", direction="up") {
  files <- list.files(path)
  for (file in files) {
    bed <- read.delim(paste(path, file, sep="/"), header=FALSE)
    bed_out <- shiftBedPositions(bed, shift=shift, pos=pos, direction=direction)
    write.table(bed_out, file=paste(path, paste(file, pos, direction, shift, sep="_"), sep="/"),
                quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
  }
}

randomizeBed <- function(data, name, write=TRUE) {
  #data <- read.delim(bed, header=F, colClasses=bedClasses)

  #for each bed record
  #choose chr based with probability based on chromosome size fraction
  #choose random position within chromosome
  #define new window based on size of old window
  #write out
  chr_lengths <- read.delim("/seq/lib/mm9_chromosome_sizes", header=F, colClasses=c("character", "numeric", "numeric"))
  chr_lengths <- chr_lengths[-nrow(chr_lengths),]
  #print(chr_lengths)
  out <- foreach (record=isplitRows(data, chunks=6), .combine="rbind") %dopar% {
    subout <- foreach (i=isplitRows(record, chunkSize=1), .combine="rbind") %do% {
      chr <- sample(chr_lengths[,1], 1, prob=chr_lengths[,3])
      start <- sample(1:chr_lengths[grep(chr, chr_lengths[,1]), 2], 1)
      span <- i[,3] - i[,2]    
      end <- start + span
      return(c(chr, start, end, i[,4], i[,5], i[,6]))
    }  
    return(subout)
  }
  if (write) {
    output_dir <- paste(bed, "random", sep="_")
    dir.create(output_dir, showWarnings=FALSE)
    write.table(out, file=paste(output_dir, name, sep="/"), quote=F, sep="\t", row.names=F, col.names=F)
  } else {
    return(out)
  }  
}

randomizeBed.batch <- function(bed, N=100) {
  pb <- txtProgressBar(min = 0, max = length(seq), style=3)
  for (i in 1:N) {
    setTxtProgressBar(pb, i) 
    randomizeBed(bed, i)
  }
}

#take data.frame with sequence names and sequence, and generate
#fasta file
formatFasta <- function(data, fname=NULL) {
  names <- as.character(data[,1])
  names <- lapply(names, function(x) gsub(" ", "", x))
  seq <- as.character(data[,2])
  fc <- file(fname, 'w')
  for (i in 1:nrow(data)) {
    name.out <- paste(">", names[i], "\n", sep="")
    cat(name.out, file=fc)
    seq.out <- paste(c(seq[i], "\n"), sep="")
    cat(seq.out, file=fc)
  }
  close(fc)
}
# For each element in BEDA find distance each element in BEDB
bedDistances <- function(a, b, chrs=NULL) {
  asplit <- split(a, a[,1])

  bsplit <- split(b, b[,1])
  if (is.null(chrs)) {
    chrs <- intersect(names(asplit), names(bsplit))
  }  
  a <- do.call("rbind", asplit[chrs])
  b <- do.call("rbind", bsplit[chrs])
#  return(b)
  result4 <- foreach (chr=chrs, .inorder=TRUE) %dopar% {
    print(chr)
#    return(asplit[[chr]])
    result3 <- foreach (arow=isplitRows(asplit[[chr]], chunkSize=1), .combine="rbind") %do% {
 #     return(arow)
      result2 <- foreach (brow=isplitRows(bsplit[[chr]], chunkSize=1), .combine="c") %do% {
#        return(1)
        a_start <- arow[2]
        a_end <- arow[3]
        b_start <- brow[2]
        b_end <- brow[3]
        result <- 0
#        print(a_end)
#        print(b_start)
        if (a_end < b_start) {
          result <- b_start - a_end
        } else if (b_end < a_start) {
          result <- a_start - b_end
        }
        return(result) 
        }
      #result2 <- unlist(result2)  
      return(result2)
    }
    rownames(result3) <- asplit[[chr]][,4]
    colnames(result3) <- bsplit[[chr]][,4]

    return(as.matrix(result3))
  }
  #return(result4)
  names(result4) <- chrs
  return(result4)
}

getSeq.bed <- function(bed, extend=0) {
  if (extend > 0) {
    bed <- trimBed(bed, -1*extend)
    bed <- trimBed(bed, extend, "down")
  }
  seq <- getSeq(Mmusculus, bed[,1], bed[,2], bed[,3], strand=bed[,6], as.character=TRUE)
  names(seq) <- bed[,4]
  return(seq)
}

getSeq.masked <- function(bed) {
  bed <- na.omit(bed)
  chrs <- unique(bed[,1])
  print(chrs)
  bed.split <- split(bed, bed[,1])
  names(bed.split) <- chrs
  total.seq <- foreach(chr=chrs, .combine="c") %dopar%  {
    chr <- as.character(chr)
    print(chr)
    chr.seq <- Mmusculus[[chr]]
    active(masks(chr.seq))['RM'] <- TRUE
    bed.curr <- bed.split[[chr]]
    bed.seq <- apply(bed.curr, 1, function(line) {
      tryCatch(return(getSeq.single(chr.seq, line)), error=function(e) return(NA))
      #if (class(result) == "try-error") {
      #  return(NA)
      #} else {
      #  return(result)
      #}
    })
    names(bed.seq) <- bed.curr[,4]
    return(bed.seq)
  }
  return(total.seq)
}

getSeq.single <- function(seq, line) {
  line.seq <- subseq(seq, start=as.numeric(line[2]), end=as.numeric(line[3]))
  line.seq <- as.character(line.seq)
  line.seq <- str_replace_all(line.seq, "#", "N")
  return(line.seq)
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
  return(bed[bed[,3] > bed[,2],])
  #return(bed)
}

splitAndSave <- function(data, chunks, fname) {
  count <- 0
  a <- foreach(c=isplitRows(data, chunks=chunks)) %do% {
    count <- count + 1
    write.table(c, file=paste(fname, count, sep="_"), quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
  }
}

shiftGetSeq <- function(bed, extend, shift_seq = "TAATTA") {
    seq <- getSeq(Mmusculus, bed[,1], bed[,2], bed[,3], strand=bed[,6], as.character=TRUE)
#  return(seq)
  #Find position of core sequence with shift_seq
  seq.shift <- lapply(seq, function(x) matchPattern(shift_seq, x, max.mismatch=2))
  seq.shift.mat <- lapply(seq.shift, as.matrix)
#  return(seq.shift.mat)
#  print("last instance")
  #If more than one take last
  for (i in 1:length(seq.shift.mat)) {
    if (nrow(seq.shift.mat[[i]]) > 1) {
      seq.shift.mat[[i]] <- seq.shift.mat[[i]][nrow(seq.shift.mat[[i]]), ]
    }
  }
  
  seq.shift.mat <- do.call("rbind", seq.shift.mat)
  print("shift")
  #Shift bed by offset
  for (i in 1:nrow(seq.shift.mat)) {
    if (bed[i,6] == "+")  {
      bed[i, 3] <- bed[i, 3] + seq.shift.mat[i, 1]
      bed[i,2] <- bed[i, 2] + seq.shift.mat[i, 1]
    } else {
      bed[i, 3] <- bed[i, 3] - seq.shift.mat[i, 1]
      bed[i,2] <- bed[i, 2] - seq.shift.mat[i, 1]
    }  
  }
  #return(bed)
  print("extend")
  bed <- trimBed(bed, -1*extend)
  bed <- trimBed(bed, extend, "down")
  #return(bed)
  print("get seq 2")
  seq <- getSeq(Mmusculus, bed[,1], bed[,2], bed[,3], strand=bed[,6], as.character=TRUE)
  return(seq)
}

## Convert isodistance BED file to matrix of coded nucleotides (observation x position)
## Input: bed (BED file)
## Input: extend (symmetrically extend BED records)
## Input: shift_seq (find first occurrence of given sequence and centralize it)
## Return: observation x position matrix with nucleotides (A = 1, C = 2, G = 3, T = 4)
bed2FaMatrix <- function(bed, extend=0, shift_seq = "TAATTA") {
  seq <- list()
  if (!is.null(shift_seq)) {
    seq <- shiftGetSeq(bed, extend, shift_seq)
  } else {
    seq <- getSeq.bed(bed, extend)
  }  
  seq <- str_split(seq, "")
  mat <- do.call("rbind", seq)
  mat <- as.data.frame(mat[,2:ncol(mat)])
  rownames(mat) <- as.character(bed[,4])
  mat <- apply(mat, 2, codeDNA)
  return(mat)
}

codeDNA <- function(data) {
  Aind <- data == "A"
  Cind <- data == "C"
  Gind <- data == "G"
  Tind <- data == "T"
  data[Aind] = 1
  data[Cind] = 2
  data[Gind] = 3 
  data[Tind] = 4
  return(as.numeric(data))
}


mononuc <- expand.grid(c("A", "C", "G", "T"))
mononuc <- apply(mononuc, 1, function(x) paste(as.character(x), collapse=""))
mononuc <- cbind(mononuc, 1:length(mononuc))

dinuc <- expand.grid(c("A", "C", "G", "T"), c("A", "C", "G", "T"))
dinuc <- apply(dinuc, 1, function(x) paste(as.character(x), collapse=""))
dinuc <- cbind(dinuc, 1:length(dinuc))

trinuc <- expand.grid(c("A", "C", "G", "T"), c("A", "C", "G", "T"), c("A", "C", "G", "T"))
trinuc <- apply(trinuc, 1, function(x) paste(as.character(x), collapse=""))
trinuc <- cbind(trinuc, 1:length(trinuc))

nuc_sets = list("mono"=mononuc, "di"=dinuc, "tri"=trinuc)

expectedFrequencies <- function(monofreq, nuc_set) {
  nuc_names <- lapply(str_split(nuc_set[,1], ""), function(x) x[2:length(x)])
  freq <- lapply(nuc_names, function(x) {
    prod(monofreq[match(x, names(monofreq))])
  })
  freq <- do.call("c", freq)
  names(freq) <- nuc_set[,1]
  return(freq)
}

## Compute frequencies of nucleotide combinations by position
## Input: seq_list (list of sequences)
## Input: nuc_set (matrix of all possible nucleotide combinations for given pattern length)
## Return: Matrix of frequeicnes 
computeFrequencies <- function(seq_list, nuc_set=trinuc, norm=TRUE, fname=NULL) {
  
  mindex <- foreach (i=1:nrow(nuc_set)) %dopar% {
    print(i)
    ms <- lapply(seq_list, function(x) {
      m <- matchPattern(nuc_set[i,1], x)
      m <- as.matrix(m)[,1]
    })
#    m <- matchPattern(nuc_set[i,1], seq)
#    m <- as.matrix(m)[,1]
    ms <- do.call("c", ms)
    ms <- factor(ms, levels=1:nchar(seq_list[[1]]))
    b <- gc()
    return(ms)
  }
 
  names(mindex) <- nuc_set[,1]
  mindex <- lapply(mindex, table)
  mindex <- do.call("rbind", mindex)
  mindex <- apply(mindex, 2, function(x) x/sum(x))
  mindex <- t(na.omit(t(mindex)))

  if (norm) {
    print("Normalizing...")
    mono <- na.omit(computeFrequencies.count(seq_list, mononuc))
    mono <- apply(mono, 2, mean)
    expfreq <- expectedFrequencies(mono, nuc_set)
    mindex <- mindex / expfreq
  }
  #return(mindex)
  if (!is.null(fname)) {
    write.table(mindex, file=fname, quote=F, sep="\t")
  }
  return(mindex)
}

computeFrequencies.set <- function(seq_list, norm=TRUE, fname=NULL) {
  mat <- lapply(1:length(nuc_sets), function(x) {
    a <- gc()
    name <- names(nuc_sets)[x]
    nuc_set <- nuc_sets[[x]]
    print(name)
    fname <- paste(nuc_seq_path, paste(fname, name, sep="_"), sep="/")
    print(paste("Saving to ", fname, sep=""))
    if (file.exists(fname)) {
      print("File exists. Skipping.")
      return(0)
    }  
    if (name == "mono") norm <- FALSE
    return(computeFrequencies(seq_list, nuc_set=nuc_set, norm=norm, fname=fname))
  })
  names(mat) <- names(nuc_sets)
  return(mat)
}

computeFrequencies.count <- function(seq_list, nuc_set=trinuc) {
  mindex <- foreach (i=1:nrow(nuc_set), .combine="cbind") %dopar% {
    print(i)
    ms <- lapply(seq_list, function(x) {
      m <- countPattern(as.character(nuc_set[i,1]), x)
      #m <- as.matrix(m)[,1]
    })
    ms <- do.call("c", ms)
    #ms <- factor(ms, levels=1:nchar(seq_list[[1]]))
    return(ms)
  }
  colnames(mindex) <- nuc_set[,1]
  rownames(mindex) <- names(seq_list)
  mindex <- t(apply(mindex, 1, function(x) x/sum(x)))
  #mindex <- apply(mindex, 2, mean)
 # mindex <- apply(mindex, 2, function(x) x - mean(x))
  return(mindex)
}

countPattern.bin <- function(chr, pattern, res, wig_fc) {
  start <- 1
  end <- 1 + res
  chr_len <- length(chr)
  
  while (end <= chr_len) {
    
  }
}


countPattern.genome <- function(pattern, res, wig) {
  #params <- new("BSParams", X=Mmusculus, FUN=matchPattern)
  #bsapply(params, pattern=pattern)
  fc <- file(wig, 'w')
  chrs <- seqnames(Mmusculus)
  lengths <- seqlengths(Mmusculus)
 
  for (i in 1:length(chrs)) {
    print(chrs[i])
    
    start <- 1
    end <- res
    out_vector <- vector("numeric", length=round(lengths[i] / res))
   
    for (j in 1:length(out_vector)) {
#    while (end <= lengths[i]) {
      if (!(round(Mestart, -1) %% round(lengths[i] * .001, -1))) print(start)
      seq <-  getSeq(Mmusculus, chrs[i], start, end, as.character=TRUE)
#      print(seq)
      out_vector[j] <- countPattern(pattern, seq)
#      cat(paste(val, "\n", sep=""), file=fc)
      start <- start + res
      end <- end + res
    }
    name.out <- paste("fixedStep",
                      paste("chr", chrs[i], sep="="),
                      paste("step", res, sep="="),
                      paste("span", res, sep="="),
                      "\n", sep=" ")
    cat(name.out, file=fc)
    cat(out_vector, file=fc, sep="\n")
  }
  close(fc)

  
}
prop_coeff <- data.frame(pattern=dinuc[,1],
                         coeff=c(-17.3, -6.7, -14.3, -16.9,
                                 -8.6, -12.8, -11.2, -8.6,
                                 -15.1, -11.7, -12.8, -6.7,
                                 -11.1, -15.1, -8.6, -17.3))

slide_coeff <- data.frame(pattern=dinuc[,1],
                          coeff=c(-0.03, -0.13, 0.47, -0.37,
                                  1.46, 0.6, 0.63, 0.47,
                                  -0.07, 0.29, 0.6, -0.13,
                                  0.74, -0.07, 1.46, -0.03))
nuc_test <- data.frame(c("AAAA", "AAAT", "AAGT", "AATA", "AATT", "AGAA", "ATAA", "ATAT", "ATTA", "GAAA", "TATA"), 1:11)

# For each sequence, compute average structural features [propeller, slide]
# using coefficients:
#     [propeller] http://srs6.bionet.nsc.ru/srs6bin/cgi-bin/wgetz?-id+ZGe1gYjVS+-e+[PROPERTY:'P0000030']
#     [slide] http://srs6.bionet.nsc.ru/srs6bin/cgi-bin/wgetz?-id+ZGe1gYjVS+-e+[PROPERTY:'P0000029']
# and dinucleotide frequencies
# Return: matrix of scores
computeStruct <- function(seq_list) {
  ## Compute dinucleotide frequencies
  di <- computeFrequencies.count(seq_list, dinuc)
  prop_scores <- di %*% prop_coeff[,2]
  slide_scores <- di %*% slide_coeff[,2]
  return(data.frame(prop=prop_scores, slide=slide_scores, row.names=names(seq_list)))
}

computeNucFeatures <- function(seq_list) {
  nuc_freq <- computeFrequencies.count(seq_list, nuc_test)
  nuc_struct <- computeStruct(seq_list)
  gc <- computeFrequencies.count(seq_list, mononuc)
  gc <- apply(gc, 1, function(x) sum(x[2:3]))
  out <- cbind(gc, nuc_freq, nuc_struct)
  rownames(out) <- names(seq_list)
  return(out)
}

computeNucPeriod <- function(seq_list, nuc_set=trinuc, max_dist=2000) {
  mindex <- foreach (i=1:nrow(nuc_set)) %dopar% {
    print(i)
    dists <- lapply(seq_list, function(x) {
      
      m <- matchPattern(nuc_set[i,1], x)
      m <- as.matrix(m)[,1]
      d <- dist(m)
      d <- c(d)
      return(d)
    })
    dists <- do.call("c", dists)
    dists <- dists[!dists>max_dist]
#    return(dists)
    h <- list(breaks=0, density=0)
    if (length(dists) > 0) {
      h <- hist(dists, density=TRUE, breaks=1:max_dist, plot=FALSE)
    }  
    out <- cbind(h$breaks, h$density)
    return(out)
    
  }
  names(mindex) <- nuc_set[,1]
  return(mindex)
}
  
plotNucFreq <- function(data, ...) {
  nr <- nrow(data)
  nc <- ncol(data)
  names <- rownames(data)
  cols <- rainbow(nr, s=.75, v=.75)

  #x11()
  par(mfrow=c(nr/4, 4), mar=c(2,2,2,2), oma=c(3,3,2,2))
  a <- sapply(1:nr, function(x) {
    plot(1:nc - (nc/2) , data[x,], type="l", col=cols[x], main=names[x], xlab="", ylab="", ...)
    abline(v=0, lty=2)
  })
  mtext("Basepairs from dyad", 1, outer=TRUE, line=1)
  mtext("Normalized frequency", 2, outer=TRUE, line=1)
}

readSplitGetNames <- function(path, filter=NULL) {
  files <- list.files(path)
  if (!is.null(filter)) files <- files[grep(filter, files)]
  out <- foreach(file=files) %do% {
    data <- read.delim(paste(path, file, sep="/"), header=FALSE)
    return(data[,4])
  }
  names(out) <- files
  return(out)
}

sortByName <- function(data, ind_pos=2) {
  data_names <- names(data)
  data_ind <- unlist(lapply(data_names, function(x) as.numeric(str_split(x, "_")[[1]][ind_pos])))
  #return(data_ind)
  return(data[order(data_ind)])  
}

classifyByQuantiles <- function(vals, probs) {
  N <- seq(1:length(probs)) - 1
  qs <- quantile(vals, probs)
  cl <- vector("numeric", length=length(vals))
  for (i in 1:(length(probs)-1)) {
    cl[vals >= qs[i] & vals < qs[i+1]] <- N[i]
  }
  return(cl)
  }
getSeqByStrand <- function(bed) {
  seq <- getSeq(Mmusculus, bed[,1], bed[,2], bed[,3])
  ind <- bed[,6] == "-"
  seq <- DNAStringSet(seq)
  seq[ind] <- reverseComplement(seq[ind])
  return(seq)
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


# Calculate the occurence frequency of each base 
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

# Input exonStart/End file
## Format: chr strand starts stops name2
# Extract sequences and paste together for given gene
# Return list of sequences, with name2 for names
stitchExons <- function(genes) {
  genes_seq <- foreach(gene=isplitRows(genes, chunkSize=1)) %do% {
    #print(gene)
    #return(gene[3])
    starts <- as.numeric(unlist(strsplit(as.character(gene[3]), ",")))
    #print(starts)
    ends <- as.numeric(unlist(strsplit(as.character(gene[4]), ",")))
    exon_count <- length(starts)
    #print(exon_count)
    gene_seq <- foreach (exon_ind=exon_count, .combine="c") %do% {
      return(getSeq(Mmusculus, as.character(gene[1]), starts[exon_ind], ends[exon_ind]))
    }
    gene_seq <- paste(gene_seq, collapse="")
    return(gene_seq)
  }
  names(genes_seq) <- genes[,5]
  return(genes_seq)
}

# To input codon usage data from http://www.kazusa.or.jp/codon
# Read in codon usage file:
# Header line: starts with >
# Codon count line:
# Realign DNAStringSet by strand
CODON_USAGE <- c("CGA", "CGC", "CGG", "CGU", "AGA", "AGG", "CUA", "CUC", "CUG", "CUU", "UUA", "UUG", "UCA", "UCC", "UCG", "UCU", "AGC", "AGU", "ACA", "ACC", "ACG", "ACU", "CCA", "CCC", "CCG", "CCU", "GCA", "GCC", "GCG", "GCU", "GGA", "GGC", "GGG", "GGU", "GUA", "GUC", "GUG", "GUU", "AAA", "AAG", "AAC", "AAU", "CAA", "CAG", "CAC", "CAU", "GAA", "GAG", "GAC", "GAU", "UAC", "UAU", "UGC", "UGU", "UUC", "UUU", "AUA", "AUC", "AUU", "AUG", "UGG", "UAA", "UAG", "UGA")
readUsage <- function(fname) {
  data <- scan(fname, what=character(), sep="\n")

  ## Extract and split headers for id
  header_ind <- as.logical(sapply(data, function(x) grep(">", x)))
  #headers <- na.omit(data[header_ind])
  #head_split <- lapply(headers, strsplit, "product=\"")
  #head_split <- lapply(head_split, function(x) strsplit(x[[1]][[2]], "\\\"/protein_id="))
  #id <- unlist(lapply(head_split, function(x) x[[1]][1]))

  ## Extract and split codon usage counts into matrix
  usage <- na.omit(data[is.na(header_ind)])
  usage_split <- lapply(usage, strsplit, " ")
  usage_split <- lapply(usage_split, function(x) as.numeric(do.call("c", x)))
  usage_matrix <- do.call("rbind", usage_split)
  colnames(usage_matrix) <- CODON_USAGE
  #return(usage_matrix)

  usage_freq_matrix <- apply(usage_matrix, 1, function(x) x / sum(x))
  return(usage_freq_matrix)  
}

alignByStrand <- function(dna_set, strand) {
  minus_ind <- strand == "-"
  dna_set[minus_ind] <- reverseComplement(dna_set[minus_ind])
  atg_pattern <- vmatchPattern("ATG", dna_set)
  atg_pattern_starts <- startIndex(atg_pattern) 
  atg_pattern_lengths <- unlist(lapply(atg_pattern_starts, length))
  atg_present_ind <- atg_pattern_lengths > 0
  atg_pattern_present <- atg_pattern_starts[atg_present_ind]
  dna_string_atg <- dna_set[atg_present_ind]
  first_atg <- unlist(lapply(atg_pattern_present, function(x) x[1]))
  dna_string_trim <- lapply(1:length(dna_string_atg), function(x) {
    substr(dna_string_atg[[x]], first_atg[[x]], length(dna_string_atg[[x]]))})
  dna_string_trim_char <- unlist(lapply(dna_string_trim, toString))
  dna_string_trim_set <- DNAStringSet(dna_string_trim_char)
  return(dna_string_trim_set)
}

# Take DNAStringSet and matched strand
# Translate to protein, reverse complementing minus strand sequences first
translateWithStrand <- function(dna_set, strand) {
  dna_set_aligned <- alignByStrand(dna_set)
  prot_string <- translate(dna_set_aligned)
  return(prot_string)
}

# Take AAStringSet and compute AA usage frequencies
protAlphabetFrequency <- function(prot) {
  letter_counts <- alphabetFrequency(prot)
  letter_counts <- letter_counts[,66:91]
  colnames(letter_counts) <- LETTERS
  aa_ind <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
  letter_counts <- letter_counts[,aa_ind]
  letter_freqs <- apply(letter_counts, 1, function(obs) obs / sum(obs))
  return(letter_freqs)
}

# Take DNAStringSet and compute codon usage
CODONS <- names(GENETIC_CODE)
codonUsage <- function(dna_set) {
  codon_counts <- vector("numeric", length=64)
 # return(codon_counts)
  names(codon_counts) <- codons
 # return(codon_counts)
  for (dna_ind in 1:length(dna_set)) {
    dna <- dna_set[[dna_ind]]
    ind <- seq(1, length(dna), 3)
    ind <- ind[-length(ind)]
    for (i in ind) {
      codon <- toString(subseq(dna, start=i, end=i+2))
      codon_counts[codon] <- codon_counts[codon] + 1
    }
  }
  #return(codon_counts)
  codon_freqs <- codon_counts / sum(codon_counts)
  return(codon_freqs)
}

codonUsageByPosition <- function(dna_set) {
  codon_pos_counts <- matrix(0, ncol=64, nrow=150, dimnames=list(1:150, CODONS))
  for (dna_ind in 1:length(dna_set)) {
    dna <- dna_set[[dna_ind]]
    ind <- seq(1, by=3, length.out=nrow(codon_pos_counts))
    if (length(dna) < ind[length(ind)]) next
    #print(ind)
    ind <- ind[-length(ind)]
    for (i in ind) {
      codon <- toString(subseq(dna, start=i, end=i+2))
      #print(codon)
      matrix_ind <- (i-1)/3 + 1
      codon_pos_counts[matrix_ind,codon] <- codon_pos_counts[matrix_ind, codon] + 1
    }
  }
  return(codon_pos_counts)
}

removeClusteredGenes <- function(bed, N) {
  id <- bed[,4]
  id_s <- sapply(id, function(x) str_replace(x, "[0-9]+$", ""))
  id_table <- table(id_s)
  #return(id_table)
  id_table <- names(id_table)[id_table <= N]
  id_index <- id_s %in% id_table
  bed <- bed[id_index,]
  id_mir <- grep("Mir", bed[,4])
  return(bed[-id_mir,])
}

get
