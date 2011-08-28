setClass(Class = 'AnnoSet',
         representation = representation(
           name='character',
           set_chr='character',
           set_pos='numeric',
           set_names='character',
           set_group='numeric',
           set_strand='character',
           #anno_data='big.matrix',
           chr_names='character',
           chr_lengths='numeric',
           bin_size='numeric',
           genome_name='character',
           coupling_vectors='list',         # CouplingVector
           read_data='list',                 # ReadSet
           pattern_data='list'              # PatternSet
           ),
         prototype = prototype(
           #coupling_vectors=NULL,
           #seq_data=NULL,
           #pattern_data=NULL
           ),
         validity = function(object){
           if(FALSE) stop("")
           return(TRUE)
         }
         )

setClass(Class = "ReadSet",
         representation = representation(
           name='character',
           raw='numeric',
           #rpm='numeric',
           ams='list'                        # NormData
           ),
         prototype = prototype()
         )

setClass(Class = 'NormalizedData',
         representation = representation(
           coupling_vector_type='character',
           data='numeric',
           cali_chr='character'
           ),
         prototype = prototype()
         )

setClass(Class = "PatternSet",
         representation = representation(
           pattern='character',
           pattern_chr='character',
           pattern_pos='numeric',
           coupling_vectors='list'   #CouplingVector
           ),
         prototype = prototype()
         )

setClass(Class = "CouplingVector",
         representation = representation(
           setCoup='numeric',
           fragmentLength='numeric',
           distFunction='character'
           ),
         prototype = prototype()
         )

setClass(Class = '')


createAnnoSet <- function(annotation) {
  a <- AnnoSet.initialize(annotation, 100, BSgenome="BSgenome.Mmusculus.UCSC.mm9")
  a <- AnnoSet.getPositions(a, pattern="CG")
  a <- AnnoSet.couplingVector(a)
  return(a)
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

addDataSet <- function(AnnoSet, dataset_name, vals) {
  DataSetObj <- new('DataSet', sample = dataset_name, raw = vals)
  AnnoSet@seq_data <- c(AnnoSet@seq_data, DataSetObj)
  return(AnnoSet)
}


