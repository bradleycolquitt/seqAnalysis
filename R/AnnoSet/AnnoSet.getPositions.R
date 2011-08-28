#####
#Function that identifies genomic positions of a given sequence pattern within the specified genome.
#####

AnnoSet.getPositions <-function(data=NULL, pattern=NULL){
	
	#if(class(data)!="AnnoSet") stop("Must specify a AnnoSet object.")
	if(is.null(pattern)){stop("Must specify a sequence pattern like CG.")}	
	
	currentchr=NULL
	currentstart=NULL
	
	organism=ls(paste("package:", data@genome_name,sep=""))
	genomedata=get(organism)
        chromosomes = data@chr_names
        #chromosomes=chr_names(data)
	#anno_data <- data@anno_data
	#starts <- split(data$set_pos, data$set_chr)
        #ends <- starts + data$bin_size - 1
	##Function which creates position files by accessing Biostrings
	#Function which creates the matrix  
	total=length(chromosomes)
        pb <- txtProgressBar(min = 0, max = total, style = 3)
	for(chromosome in chromosomes){
	        #{subject<- Views(genomedata[[chromosome]], start = starts, end = ends)}
          {subject <- genomedata[[chromosome]]}
                setTxtProgressBar(pb, which(chromosomes==chromosome))
		plus_matches<-matchPattern(pattern,subject)
		start=start(plus_matches)   
		currentchr=c(currentchr,rep(chromosome,length(start)))
		currentstart=c(currentstart,start) 

		patternDNA<-DNAString(pattern)
		rcpattern<-reverseComplement(patternDNA)
		if(patternDNA!=rcpattern){
			minus_matches<-matchPattern(rcpattern,subject)
			start=start(minus_matches) 
			currentchr=c(currentchr,rep(chromosome,length(start)))
			currentstart=c(currentstart,start) 
		}			
	}

        AnnoSetObj <- .addPattern(data, pattern, currentchr, currentstart)
        cat("\n")
        return(AnnoSetObj)
        
        
# 	AnnoSetObj = new('AnnoSet', seq_pattern=as.character(pattern), pattern_chr=currentchr, pattern_pos=currentstart, number_pattern=#length(currentstart), genome_chr=genome_chr(data), genome_pos=genome_pos(data), genome_raw=genome_raw(data), extend=extend(data), bin_s#ize=bin_size(data), sample_name=sample_name(data), genome_name=genome_name(data), regions_chr=regions_chr(data), regions_start=regions_#start(data), regions_stop=regions_stop(data), regions_strand=regions_strand(data), number_regions=number_regions(data), chr_names=chr_n#ames(data), chr_lengths=chr_lengths(data))
#	cat("\n")

#	return(AnnoSetObj)	
}

.addPattern <- function(data, pattern, currentchr, currentstart) {
  PatternSetObj <- new('PatternSet', pattern=pattern, pattern_chr=currentchr, pattern_pos=currentstart)
  if(length(data@pattern_data) == 0) {
    data@pattern_data <- list(PatternSetObj)
  }
  else {
    data@pattern_data <- c(data@pattern_data, PatternSetObj)
  }  
    return(data)
}
