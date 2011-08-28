###############
##Function takes a MEDIPSset and weights all signals by the estimated mean signal of its coupling factor,
##transforms the normalized signals into a reads/million scale and finally transforms the data range into the interval of [0:1000].
###############

##Perform full normalization on a given MEDIPS.SET.
AnnoSetSimp.normalize <- function(data=NULL, reads=NULL, transf=FALSE){
	
	if(class(data)!="AnnoSetSimp") stop("Must specify a AnnoSetSimp object.")
        read_name <- strsplit(reads, ".bed")[[1]]
        cali_name <- paste(read_name, data@bin_size, data@pattern_data@pattern,
                           data@coupling_vector@distFunction, sep="_") 
        calibration_info <- loadCalibrationInfo(cali_name)
        curr_read_set <- .extract(data@read_data, reads)
        signal <- curr_read_set@raw
        coupling <- data@coupling_vector@setCoup
        intercept <-calibration_info[1]
        slope <- calibration_info[2]           
	##Weight signals by linear regression obtained parameters
	####################
	cat("Weight raw signals by the estimated parameters...\n")
	estimated_meanSignal=(coupling-intercept)/slope	
	if(length(estimated_meanSignal[estimated_meanSignal<=0])!=0){
		signal[estimated_meanSignal>0]=signal[estimated_meanSignal>0]/estimated_meanSignal[estimated_meanSignal>0]
		minWeight=min(estimated_meanSignal[estimated_meanSignal>0])
		signal[estimated_meanSignal<=0]=signal[estimated_meanSignal<=0]/minWeight
	}
	else{signal=signal/estimated_meanSignal}		
	
	##Transform weighted signals into reads per million
	####################
	cat("Transform weighted signals into reads per million (rpm)...\n")
	signal <- signal/curr_read_set@number_regions/1000000
        if (transf) signal <- MEDIPS.transform(signal)
        .updateReadData <- function(data) {
          ReadDataObj <- new('ReadSet', name=curr_read_set@name, raw=curr_read_set@raw,
                             rms=signal, number_regions=curr_read_set@number_regions)
          data@read_data <- .replace(data@read_data, reads, ReadDataObj)
          return(data)
        }
        
        return(.updateReadData(data))

	AnnoSetObj = new('AnnoSet', genome_norm=signal, cali_chr=cali_chr(data), calcurve_mean_signals=calcurve_mean_signals(data), calcurve_mean_coupling=calcurve_mean_coupling(data), calcurve_var=calcurve_var(data), intercept=intercept(data), slope=slope(data), genome_CF=genome_CF(data), fragmentLength=fragmentLength(data), distFunction=distFunction(data), distFile=distFile(data), seq_pattern=seq_pattern(data), pattern_chr=pattern_chr(data), pattern_pos=pattern_pos(data), number_pattern=number_pattern(data), genome_chr=genome_chr(data), genome_pos=genome_pos(data), genome_raw=genome_raw(data), extend=extend(data), bin_size=bin_size(data), sample_name=sample_name(data), genome_name=genome_name(data), regions_chr=regions_chr(data), regions_start=regions_start(data), regions_stop=regions_stop(data), regions_strand=regions_strand(data), number_regions=number_regions(data), chr_names=chr_names(data), chr_lengths=chr_lengths(data))
}

##Only transform given data into log scale and shift into interval
##################################################################
MEDIPS.transform <- function(data=NULL){
	##Log2 of signals except signal=0
	####################
	data[!is.na(data)] = log2(data[!is.na(data)])
				
	##Shift signals into positive value range
	####################
	minsignal=min(data[data!=-Inf & !is.na(data)])
	data[data!=-Inf & !is.na(data)]=data[data!=-Inf & !is.na(data)]+abs(minsignal)
	
	##Transform values into the interval [0:1000]
	####################
	maxsignal=max(data[data!=Inf & !is.na(data)])
	data[!is.na(data)]=(data[!is.na(data)]/maxsignal)*1000

	##Eliminate -Inf --> 0
	#######################
	data[data==-Inf]=0

	gc()
	return(data)
}

