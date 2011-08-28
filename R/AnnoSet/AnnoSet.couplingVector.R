##################
##The function calculates the pattern influence factor for all given pattern positions and all genomic bins.
##It returns a vector that lists for each genomic bin the total pattern coupling factor.
##################

library(MEDIPS)
#dyn.load("/usr/lib64/R/lib/coupling.so")


AnnoSet.couplingVector <- function(data=NULL, pattern="CG", distFile="empty", fragmentLength=700, func="count"){
	
#	if(class(data)!="MEDIPSset") stop("Must specify a MEDIPSset object.")	
	if (func!="linear" && func!="exp" && func!="log" && func!="count" && func!="custom"){stop("func has to be linear, exp, log, count, or custom.")}
	if(func=="custom" && distFile=="empty"){stop("For the custom distance function, a distance file has to be stated!")}
	if(func=="custom" && distFile!="empty"){
		print(paste("Reading custom distance file ", distFile, sep=" "), quote=F)	
		distFileData=read.table(distFile, sep='\t', header=F, row.names=NULL, colClasses=c("numeric", "numeric"))
	}
        
        set_chr = data@set_chr
        set_pos = data@set_pos
        current_pattern <- .extract(data@pattern_data, pattern)
        pattern_chr = current_pattern@pattern_chr
        pattern_pos = current_pattern@pattern_pos
        chromosomes = data@chr_names
        chr_lengths = data@chr_lengths
        bin_size = data@bin_size
        
#	genome_chr=genome_chr(data)
#	genome_pos=genome_pos(data)
#	pattern_chr=pattern_chr(data)
#	pattern_pos=pattern_pos(data)	
#	chromosomes=chr_names(data)
#	chr_lengths=chr_lengths(data)
#	bin_size=bin_size(data)
		
	##Initialization of the genome coupling vector.	
	setCoup=vector(length=length(set_chr), mode="numeric")	
	setCoup[]=0
	
	##Define coupling factor functions.
	if(func=="linear"){wfun=function(dista, FLength){return(1-dista/FLength)}}	
	if(func=="exp"){wfun=function(dista, FLength){return(1 - dista^2/(FLength)^2)}}
	if(func=="log"){wfun=function(dista, FLength){return(1 - log(1 + abs(dista)/(FLength/18), 10))}}    
	if(func=="count"){wfun=function(dista, FLength){return(1)}}
	if(func=="custom"){wfun=function(dista, FLength){if(dista>(length(distFileData[,1])-1)){return(0)}else{return(distFileData[(dista+1),2])}}}
				
	##Pre-compute coupling distance weights.
	distanceVector=NULL
	for(i in 0:(fragmentLength-1)){
		distanceVector=as.numeric(c(distanceVector, wfun(i, fragmentLength)))
	}
	##Calculate genome wide coupling factors	
	cat("Calculating coupling factors...\n")
	noInfluenceFlankWin=floor(fragmentLength/bin_size)
	total=length(chromosomes)
        pb <- txtProgressBar(min = 0, max = total, style = 3)
	for(i in 1:length(chromosomes)){
#   		cat(paste("Calculating coupling factors for", chromosomes[i],"....\n", sep=" "))	
 		setTxtProgressBar(pb, i)
		temp_genomePos=set_pos[set_chr==chromosomes[i]]
		temp_setCoup=vector(length=length(temp_genomePos), mode="numeric")
		temp_setCoup[]=0	
		temp_posPos=pattern_pos[pattern_chr==chromosomes[i]]
		temp_setCoup=.Call("coupling",temp_posPos,bin_size,as.integer(temp_genomePos),temp_setCoup,noInfluenceFlankWin,fragmentLength-1,distanceVector)
 		setCoup[set_chr==chromosomes[i]]=as.numeric(temp_setCoup)			
	}

        #return(setCoup)
        
	AnnoSetObj <- .addCouplingVector(data, pattern, setCoup, fragmentLength = fragmentLength, distFunction = func)		
        cat("\n")
        return(AnnoSetObj)
#	MEDIPSsetObj = new('MEDIPSset', genome_CF=setCoup, fragmentLength=fragmentLength, distFunction=func, distFile=distFile, seq_pattern=seq_pattern(data), pattern_chr=pattern_chr(data), pattern_pos=pattern_pos(data), number_pattern=number_pattern(data), set_chr=set_chr(data), set_pos=set_pos(data), genome_raw=genome_raw(data), extend=extend(data), bin_size=bin_size(data), sample_name=sample_name(data), genome_name=genome_name(data), regions_chr=regions_chr(data), regions_start=regions_start(data), regions_stop=regions_stop(data), regions_strand=regions_strand(data), number_regions=number_regions(data), chr_names=chr_names(data), chr_lengths=chr_lengths(data))
#	cat("\n")
	
#	return(MEDIPSsetObj)	
}


.find <- function(set, name) {
  index <- unlist(sapply(1:length(set), function(x) {
    if(set[[x]]@pattern ==  name) {
      return(x)
    }
  }
  ))                
  return(as.numeric(index))   
}

.extract <- function(set, name) {
  return(set[[.find(set, name)]])
}

.replace <- function(set, name, replacement) {
  index <- .find(set, name)
  set[index] <- replacement
  return(set)
}

.addCouplingVector <- function(data, pattern, setCoup, fragmentLength, distFunction) {
  current_pattern <- .extract(data@pattern_data, pattern)
  CouplingVectorObj <- new('CouplingVector', setCoup = setCoup, fragmentLength=fragmentLength, distFunction=distFunction)
  if(length(current_pattern@coupling_vectors) == 0) {
    current_pattern@coupling_vectors<- list(CouplingVectorObj)
  }
  else {
    current_pattern@coupling_vectors <- c(current_pattern@coupling_vectors, CouplingVectorObj)
  }
  data@pattern_data <- .replace(data@pattern_data, pattern, current_pattern)
  return(data)

}
