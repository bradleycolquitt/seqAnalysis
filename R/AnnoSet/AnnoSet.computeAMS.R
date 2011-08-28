## Take AnnoSet object
## Divide RMS by coupling vector
## Transform to [0,1000]

AnnoSet.computeAMS <- function(data, reads) {
  read.data <- .extract(data@read_data, reads)
  rms <- read.data@rms
  cf <- data@coupling_vector@setCoup
  
  ams <- rms/cf
  return(cf)
  ams <- MEDIPS.transform(ams)

  .updateReadData <- function(data) {
    ReadDataObj <- new('ReadSet', name=read.data@name, raw=read.data@raw,
                       rms=read.data@rms, ams=ams, number_regions=read.data@number_regions)
    data@read_data <- .replace(data@read_data, reads, ReadDataObj)
    return(data)
        }
        
  return(.updateReadData(data))
}
