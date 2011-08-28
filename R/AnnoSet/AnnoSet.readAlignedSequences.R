#####
#Function that reads a region file.
#Region file is a simple tab-delimited text files (chr | start | stop | strand).
#####

AnnoSet.readAlignedSequences <-
function(data=NULL, fileName=NULL, numrows=-1, extend=350){
		
	## Read region file	
	regions=NULL
	path <- "/home/user/data/working"
        file <- paste(path, fileName, sep="/")
	
	cat(paste("Reading file ", fileName, " in ", path, "...\n", sep=""))		
	
	if(!fileName%in%dir(path)){
		stop(paste("File", fileName, " not found in", path, sep =" "))			
	}
	
	regions=read.table(file, sep='\t', header=FALSE, row.names=NULL, nrows=numrows,
          colClasses=c("character", "numeric", "numeric", "character"))			
		
        chr <- regions[,1]
        start <- regions[,2]
        stop <- regions[,3]
        strand <- regions[,4]
        chromosomes <- data@chr_names
        set_chr <- data@set_chr
        set_pos <- data@set_pos
        chr_split <- split(set_chr, set_chr)
        no_chr_windows <- sapply(chr_split, length)
        supersize_chr <- cumsum(no_chr_windows)
        
        ##Distribute reads over set.	
        total = length(chromosomes)
        pb <- txtProgressBar(min = 0, max = total, style = 3)
        genomeVec_signal=vector(length=supersize_chr[length(chromosomes)], mode="numeric")
	cat("\nDistribute reads over set...\n")
        for(i in 1:length(chromosomes)){
		setTxtProgressBar(pb, i)
                genomeVec_signal[set_chr==chromosomes[i]] <- MEDIPS.distributeReads(start[chr==chromosomes[i]], stop[chr==chromosomes[i]],
                                   strand[chr==chromosomes[i]], set_pos[set_chr==chromosomes[i]], extend)
	}
	cat("\n")

	.addDataSet <- function(data, signal, name) {
           ReadSetObj <- new('ReadSet', name=name, raw=signal, norm=0, number_regions=nrow(regions))
           if(length(data@read_data) == 0) {
             data@read_data <- list(ReadSetObj)
           }
           else if (.find(data@read_data, name) == 0) {
             data@read_data <- c(data@read_data, ReadSetObj)
           }
           else {
             data@read_data <- .replace(data@read_data, name, ReadSetObj)
           }
           return(data)
        }

        return(.addDataSet(data, genomeVec_signal, fileName))   
}
