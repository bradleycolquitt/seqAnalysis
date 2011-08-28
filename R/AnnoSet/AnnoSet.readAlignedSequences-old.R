#####
#Function that reads a region file.
#Region file is a simple tab-delimited text files (chr | start | stop | strand).
#####

AnnoSet.readAlignedSequences <-
function(data=NULL, fileName=NULL, numrows=-1, extend=350){
		
	## Species dependent definitions.
	#if(is.null(BSgenome)){stop("Must specify a BSgenome library.")}
	if(is.null(file)){stop("Must specify a region file.")}
		
	## Read region file	
	regions=NULL
	path <- "/home/user/data/working"
        file <- paste(path, fileName, sep="/")
	#fileName=unlist(strsplit(file, "/"))[length(unlist(strsplit(file, "/")))]
	#path=paste(unlist(strsplit(file, "/"))[1:(length(unlist(strsplit(file, "/"))))-1],collapse="/") 
#	if(path==""){
#		path=getwd()
#	}
	
	cat(paste("Reading file ", fileName, " in ", path, "...\n", sep=""))		
	
	if(!fileName%in%dir(path)){
		stop(paste("File", fileName, " not found in", path, sep =" "))			
	}
	
	regions=read.table(file, sep='\t', header=FALSE, row.names=NULL, nrows=numrows,
          colClasses=c("character", "numeric", "numeric", "character"))			
		
	## Sort chromosomes
	#if(length(unique(regions[,1]))>1){chromosomes=mixedsort(unique(regions[,1]))}
	#if(length(unique(regions[,1]))==1){chromosomes=unique(regions[,1])}
	
	## Get chromosome lengths for all chromosomes within data set.
	#cat(paste("Loading chromosome lengths for ",BSgenome, "...\n", sep=""))		
	#dataset=get(ls(paste("package:", BSgenome, sep="")))
	#chr_lengths=as.numeric(sapply(chromosomes, function(x){as.numeric(length(dataset[[x]]))}))

        chr <- regions[,1]
        start <- regions[,2]
        stop <- regions[,3]
        strand <- regions[,4]
        chromosomes <- data@chr_names
        set_chr <- data@set_chr
        set_pos <- data@set_pos
        #no_chr_windows <- ceiling(chr_lengths/data@bin_size)
        chr_split <- split(set_chr, set_chr)
        no_chr_windows <- sapply(chr_split, length)
        supersize_chr <- cumsum(no_chr_windows)
        
	## Create the genome vector
	cat("Create the genome vector...\n")
	#setVec_chr=vector(length=supersize_chr[length(chromosomes)], mode="character")
	#setVec_pos=vector(length=supersize_chr[length(chromosomes)], mode="numeric")		
	total=length(chromosomes)
        pb <- txtProgressBar(min = 0, max = total, style = 3)
        #lengths <- list()
# 	for(i in 1:length(chromosomes)){
#		setTxtProgressBar(pb, i)
#		if(i==1){
 #                       index <- c(1:no_chr_windows[i])
#			setVec_chr[index]=chromosomes[i]
			#setVec_pos[1:no_chr_windows[i]]=seq(1, chr_lengths[i], bin_size)
#                        setVec_pos[index] <- set_pos[index]
                        #lengths <- c(lengths, length(set_chr==chromosomes[i]))
#                }
#		if(i>1){
#                        index <- c((supersize_chr[i-1]+1):(supersize_chr[i-1]+no_chr_windows[i]))
#			setVec_chr[index]=chromosomes[i]
			#genomeVec_pos[(supersize_chr[i-1]+1):(supersize_chr[i-1]+no_chr_windows[i])]=seq(1, chr_lengths[i], bin_size)
#                        setVec_pos[index] <- set_pos[index]
                        #lengths <- c(lengths, length(set_chr==chromosomes[i]))
#                      }

#        }
	#return(lengths)
	##Distribute reads over genome.	
	genomeVec_signal=vector(length=supersize_chr[length(chromosomes)], mode="numeric")
	cat("\nDistribute reads over genome...\n")
        #pb <- txtProgressBar(min = 0, max = total, style = 3)
        for(i in 1:length(chromosomes)){
		setTxtProgressBar(pb, i)
		#genomeVec_signal[genomeVec_chr==chromosomes[i]] = MEDIPS.distributeReads(start[chr==chromosomes[i]], stop[chr==chromosomes[i]], strand[chr==chromosomes[i]], genomeVec_pos[genomeVec_chr==chromosomes[i]], extend)
               # genomeVec_signal[setVec_chr==chromosomes[i]] = MEDIPS.distributeReads(start[chr==chromosomes[i]], stop[chr==chromosomes[i]], strand[chr==chromosomes[i]], setVec_pos[setVec_chr==chromosomes[i]], extend)
                genomeVec_signal[set_chr==chromosomes[i]] <- MEDIPS.distributeReads(start[chr==chromosomes[i]], stop[chr==chromosomes[i]], strand[chr==chromosomes[i]], set_pos[set_chr==chromosomes[i]], extend)
	
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

	## creating AnnoSet object
#    	AnnoSetObj = new('AnnoSet', sample_name=fileName, genome_name=BSgenome, regions_chr=regions[,1], regions_start=regions[,2], regions_stop=regions[,3], regions_strand=regions[,4], number_regions=length(regions[,1]), chr_names=chromosomes, chr_lengths=chr_lengths)
#    	return(AnnoSetObj)    
}
