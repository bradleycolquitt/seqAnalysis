### MP all MEDIPS for all annotations
### NOT FUNC
MP.all <- function(bin.size=200) {
  medips.path <- "~/data/medips"
  dir_files <- list.files(medips.path)
  select_files <- dir_files[grep(as.character(bin.size), dir_files)]
  for(file in select_files) {
    load(paste(medips.path, file, sep="/"))
    
  }
}
