#####
#Function that identifies genomic positions of a given sequence pattern within the specified genome.
#####

AnnoSetSimp.getPositions <-function(data=NULL, pattern=NULL){
	
	if(is.null(pattern)){stop("Must specify a sequence pattern like CG.")}	
	
	currentchr=NULL
	currentstart=NULL
	
	organism=ls(paste("package:", data@genome_name,sep=""))
	genomedata=get(organism)
        chromosomes = data@chr_names
        #chromosomes=chr_names(data)
	#anno_data <- data@anno_data
	starts <- split(data@set_pos, data@set_chr)
        ends <- lapply(starts, function(x) x + data@bin_size - 1)
        names(starts) <- chromosomes
        names(ends) <- chromosomes
	##Function which creates position files by accessing Biostrings
	#Function which creates the matrix  
	total=length(chromosomes)
        pb <- txtProgressBar(min = 0, max = total, style = 3)
	for(chromosome in chromosomes){
	        subject <- Views(genomedata[[chromosome]],
                                start = starts[[chromosome]],
                                end = ends[[chromosome]])
          #subject <- genomedata[[chromosome]]
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
        return(currentstart)
        AnnoSetObj <- .addPattern(data, pattern, currentchr, currentstart)
        cat("\n")
        return(AnnoSetObj)
        
}

.addPattern <- function(data, pattern, currentchr, currentstart) {
  PatternSetObj <- new('PatternSet', pattern=pattern, pattern_chr=currentchr, pattern_pos=currentstart)
  data@pattern_data <- PatternSetObj
  #if(length(data@pattern_data) == 0) {
  #  data@pattern_data <- list(PatternSetObj)
  #}
  #else {
  #  data@pattern_data <- c(data@pattern_data, PatternSetObj)
  #}  
  return(data)
}
