## Identify positions of HD associated 8-mers within a given sequence
## Input sequence
## Output matrix of relative positions and E-scores

# Libraries
library(BSgenome)
library(foreach)
library(doMC)

# Source files
source("~/src/seqAnalysis/R/seqUtil.R")

# Options
registerDoMC(cores=4)


# Matrix of 8-mer E-scores from Berger, Cell, 2008
# Uncomment to load
# 32896 rows of 8-mers
# 180 columns: 1-2 are 8-mer sequence and reverse complement, 3-180 are HD
# E-scores computed by
#   - defining 'foreground' (contains given 8-mer sequence) or 'background' (does not contain)
#   - took top half of each set and computed a modified Wilcoxon-Mann-Whitney stat (essentially comparing ranks) scaled to be invariant for sample size
#Data ranged from -0.5 to 0.5
#ESCORE_MAT <- read.delim("~/data/hughes/table_hox_enrichmentscore_102107.txt")

# Loop through target sequences, collect as list
#   Matrix <- Step along target seuqence
#     Identify 8-mer
#     return E-scores (rbind)
#   return Matrix
collectEscores <- function(seq_list) {
  n <- 1
  e_list <- foreach(seq=seq_list) %dopar% {
    escore_ind <- 0
    n <- n + 1
    print(n)
    position_matrix <- foreach(i=1:(width(seq) - 7), .inorder=TRUE, .combine="rbind") %do% {
      #print(i)
      seq8 <- subseq(seq, i, i+7)
      escore_ind <- grep(seq8, ESCORE_MAT[,1])
      if (length(escore_ind) == 0) {
        escore_ind <- grep(seq8, ESCORE_MAT[,2])
      }
      return(ESCORE_MAT[escore_ind, 3:180])
    }
    n <- n + 1
    return(position_matrix)
  }
  
  return(e_list)
}

# Identify HD with greatest signal
# Loop through target-Escore matrix
#   compute column sums
#   Return N columns with highest sums
sumScores <- function(seq_e, N=10) {
  lapply(seq_e, function(x) {
    ordered_matrix <- x[,order(-apply(x, 2, sum))]
    return(ordered_matrix[,1:10])    
  })
}

filterByEscore <- function(thresh=0.45) {
  result <- apply(ESCORE_MAT[,3:180], 2, function(col) {
    ind <- col >= thresh
    return(cbind(ESCORE_MAT[ind, 1:2], col[ind]))
  })
  return(result)
}

# Determine number of HD sites (defined as 8-mer with E>=thresh) in sequence
# Input list of sequences, E-score threshold
# Output matrix with counts of HD site occurances
### USE Jellyfish to count
count_hd_sites <- function(seq_list, thresh) {
  print("thresh")
  escore_thresh <- cbind(ESCORE_MAT[,1:2], apply(ESCORE_MAT[,3:180] >= thresh, 2, as.numeric))
  print("loop1")
  # Encode 8-mer library as DNAStringSet
  mer8 <- DNAStringSet(c(as.character(ESCORE_MAT[,1]), as.character(ESCORE_MAT[,2])))
  mer_length <- nrow(ESCORE_MAT)                     
  #mer8.rc <- reverseComplement(mer8)
  # Loop through sequence list
  result <- foreach(seq=seq_list, .combine="rbind") %dopar% {
    
    # Step through sequence, match each 8-mer to thresholded E-score matrix
    position_matrix <- foreach(i=1:(length(seq) - 7), .inorder=TRUE, .combine="rbind") %do% {
      seq8 <- subseq(seq, i, i+7)
      mindex <- vmatchPattern(seq8, mer8, algorithm="naive-exact")
      escore_ind <- which(countIndex(ind)>0)
      escore_ind <- escore_ind %% mer_length
      #escore_ind <- grep(seq8, escore_thresh[,1])
      #if (length(escore_ind) == 0) {
      #  mindex <- vmatchPattern(seq8, mer8.rc, algorithm="naive-exact")
      #  escore_ind <- which(countIndex(ind)>0)
        #escore_ind <- grep(seq8, escore_thresh[,2])
      #}
      return(escore_thresh[escore_ind, 3:180])
    }
    position_matrix <- colSums(position_matrix)
    return(position_matrix)
  }
  rownames(result)<- names(seq_list)
  return(result)
}



# Read in fasta files of kmer occurences produced by jellyfish
# Fill matrix of these counts
parse_kmer_counts <- function(dir) {
  files <- list.files(dir)
  files <- files[grep("fa$", files)]
  
  files_prefixes <- str_split(files, ".fa")
  files_prefixes <- unlist(lapply(files_prefixes, function(x) x[1]))
  #files_prefixes <- unlist(lapply(files_prefixes, function(x) paste("el", x, sep="_")))
  # Rows: sequences
  # Cols: 8mer
  count_matrix <- matrix(0, nrow=length(files_prefixes), ncol=nrow(ESCORE_MAT)*2, 
                         dimnames=list(files_prefixes, c(as.character(ESCORE_MAT[,1]), as.character(ESCORE_MAT[,2]))))

  # Fill matrix with occurence information for each 8mer
  print("Filling count matrix...")
  pb <- txtProgressBar(min = 0, max = length(files_prefixes), style = 3)
  for (i in 1:length(files_prefixes)) {
    setTxtProgressBar(pb, i)
    fasta <- read.DNAStringSet(paste(dir, files[i], sep="/"), "fasta")
    count_matrix[i,match(as.character(fasta), colnames(count_matrix))] <- names(fasta)
    }
  close(pb)
  
  count_matrix <- apply(count_matrix, 2, as.numeric)
  return(count_matrix)
}
  
# Multiply count matrix by thresholded ESCORE_MAT to remove non-signifcant 8mers
threshold_count_matrix <- function(count_matrix) {
  print("Thresholding count matrix...")
  escore_thresh <- apply(ESCORE_MAT[,3:180] >= 0.45, 2, as.numeric)
  escore_thresh <- rbind(escore_thresh, escore_thresh)
  result_matrix <- count_matrix %*% escore_thresh
  return(result_matrix)
}

# Normalize count matrix by expected occurence by binomial
norm_by_binom <- function(counts, seq_list) {
  # Determine mononucleotide frequencies
  print("Compute mononucleotide frequencies")
  mononuc_freq <- computeFrequencies.count(seq_list)
  mononuc_freq_mean <- apply(mononuc_freq, 2, mean, na.rm=TRUE)
  
  # Loop across colnames (kmers) and calculate expected frequency for each
  print("Computing expected frequencies")
  kmer_freq <- foreach(kmer=colnames(counts), .combine="c") %dopar% {
    expectedFrequencies.single(kmer, mononuc_freq_mean)
  }
  names(kmer_freq) <- colnames(counts)
  
  # Construct expectation matrix from
  # vector of sequence lengths and vector of kmer frequencies
  print("Construct expectation matrix")
  kmer_freq_seq <- width(seq_list) %*% t(kmer_freq)
  rownames(kmer_freq_seq) <- names(seq_list)
  #return(kmer_freq_seq)
  
  # Normalize counts by expectation matrix
  print("Normalize counts")
  kmer_freq_seq <- kmer_freq_seq[match(rownames(counts), rownames(kmer_freq_seq)),]
  kmer_freq_seq <- na.omit(kmer_freq_seq)
  kmer_freq_seq <- kmer_freq_seq[-grep(setdiff(rownames(kmer_freq_seq), rownames(counts)), rownames(kmer_freq_seq)),]
  #return(kmer_freq_seq)
  counts_norm <- counts / kmer_freq_seq
  
  return(counts_norm)
}