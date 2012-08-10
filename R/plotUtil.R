checkExisting <- function(input.path, output.path, select=NULL, extend=NULL) {
  input.files <- list.files(input.path)
  if (!is.null(select)) input.files <- input.files[input.files%in%select]
  output.files <- list.files(output.path)
  if (!is.null(extend)) output.files <- output.files[grep(extend, output.files)]
  files.tbp <- sapply(input.files, function(x) grep(x, output.files))
  files.tbp.index <- unlist(lapply(files.tbp, function(x) is.na(x[1])))
  files.tbp <- names(files.tbp[files.tbp.index])
  return(files.tbp)
}

profileRead <- function(data.path, fun, select, select2=NULL) {
  files <- list.files(data.path)
  #print(files)
  files <- files[grep(select, files)]
  #print(files)
  if(length(files) == 0) {
    stop("Annotation file not found")
  }
  val_files <- NULL
  ci_files_1 <- NULL
  ci_files_2 <- NULL
  #files <- files[grep(select, files)]
  if (!is.null(select2)) {
    #files <- files[grep(select2, files)]
    filter <- paste(select, select2, sep="_")
    #print(files)
    #val_files <- sort(files[grep(paste(select, "mean$", sep="_"), files)])
    #ci_files_1 <- sort(files[grep(paste(select, "bootCI_1", sep="_"), files)])
    #ci_files_2 <- sort(files[grep(paste(select, "bootCI_2", sep="_"), files)])
    #val_files <- sort(files[grep(paste(filter, "mean$", sep="_"), files)])
    #ci_files_1 <- sort(files[grep(paste(filter, "mean_bootCI_1", sep="_"), files)])
    #ci_files_2 <- sort(files[grep(paste(filter, "mean_bootCI_2", sep="_"), files)])
    val_files <- sort(files[grep(paste(filter, paste(fun, "$", sep=""), sep="_"), files)])
    ci_files_1 <- sort(files[grep(paste(filter, fun, "bootCI_1", sep="_"), files)])
    ci_files_2 <- sort(files[grep(paste(filter, fun, "bootCI_2", sep="_"), files)])
  } else {
    #val_files <- sort(files[grep(paste(select, "mean$", sep="_"), files)])
    #ci_files_1 <- sort(files[grep(paste(select, "mean_bootCI_1", sep="_"), files)])
    #ci_files_2 <- sort(files[grep(paste(select, "mean_bootCI_2", sep="_"), files)])
    val_files <- sort(files[grep(paste(select, paste(fun, "$", sep=""), sep="_"), files)])
    ci_files_1 <- sort(files[grep(paste(select, fun, "bootCI_1", sep="_"), files)])
    ci_files_2 <- sort(files[grep(paste(select, fun, "bootCI_2", sep="_"), files)])
    #val_files <- sort(files[grep("mean$", files)])
    #ci_files_1 <- sort(files[grep("mean_bootCI_1", files)])
    #ci_files_2 <- sort(files[grep("mean_bootCI_2", files)])
  }
  #print(val_files)
  #print(ci_files_1)
  #print(ci_files_2)
  val_data <- foreach(file=val_files, .combine="cbind") %do% {
    read.delim(paste(data.path, file, sep="/"), header=F)}
  ci_data_1 <- foreach(file=ci_files_1 , .combine="cbind") %do% {
    read.delim(paste(data.path, file, sep="/"), header=F)}
  ci_data_2 <- foreach(file=ci_files_2, .combine="cbind") %do% {
    read.delim(paste(data.path, file, sep="/"), header=F)}
  return(list(val_data, ci_data_1, ci_data_2))   
}

MP.getData <- function(path, filters, neg_filter=NULL, annotation, group2) {
  files <- list.files(path)
  if (!is.null(neg_filter)) {
    for (nf in neg_filter) {
      files <- files[-grep(nf, files)]
    }
  }  
  for (f in filters) {
    files <- files[grep(f, files)]
  }
  files <- sort(files)
  data <- lapply(files, function(x) {
    file.path <- paste(path, x, "profiles", sep="/")
#    print(file.path)
#    print(annotation)
    profileRead(file.path, annotation, group2)
  })
  return(data)
}

splitReform <- function(data.in) { 
foreach (i=1:4) %do% {
    data <- list(data.in[[i]])
    data <- c(data, list(data.in[[i+4]]))
    data <- c(data, list(data.in[[i+8]]))
  }
}

trimData <- function(data, range) {
  data <- lapply(data, function(x) {
    lapply(x, function(y) {
      y[range[1]:range[2]]
    })
  })
  return(data)
}
getRange <- function(data, buffer=100) {
  # print(length(data))
  ci_lower <- lapply(data, function(x) x[[2]])
  ci_upper <- lapply(data, function(x) x[[3]])
  val.max <- 0
  val.min <- 1000
  #print(ci_lower)
  #print(dim(ci_lower[[1]]))
  for(i in 1:length(ci_upper)) val.max <- max(val.max, max(ci_upper[[i]]))
  for(i in 1:length(ci_lower)) val.min <- min(val.min, min(ci_lower[[i]]))
  buffer <- (val.max - val.min) / 4
  val.max <- round(val.max + buffer,3)
  val.min <- round(val.min - buffer,3)
  val <- c(val.min, val.max)
}

findMaxMin <- function(x) {
  val.max <- 0
  val.min <- 1000
  for (i in 1:length(x)) {
    val.max <- max(val.max, x[[i]][2])
    val.min <- min(val.min, x[[i]][1])
  }
  return(c(val.min, val.max))
}

baselineNorm <- function(data) {
  starts <- unlist(lapply(data, function(x) mean(x[[1]][[1]][1:10])))
  norm_val <- mean(starts)
  ind_norm_vals <- starts - norm_val
  data <- lapply(1:length(data), function(x) {
    lapply(data[[x]], function(y) y - ind_norm_vals[[x]])
  })
  return(data)
}

computeAxis <- function(data, wsize, lab) {
  lab.l <- length(lab)
  lab.data <- NULL
  l <- length(data)
#  print(l)
  a <- 1
  if (lab.l == 1) {
    b <- l / 2 + 1
    c <- l
    lab_dist <- l / 2 * wsize / 1000
    lab.data <- list(pos=c(a, b, c), dist=lab_dist)
  } else if (lab.l == 2) {
    b <- round(l - 1, digits=-2) / 2  + 1 
    c <- b + l - (b-1) * 2
    d <- l
    lab_dist <- (b - 1) * wsize / 1000
    lab.data <- list(pos=c(a, b, c, d), dist=lab_dist)
 #   print(lab.data)
  }
}

.unwrap <- function(list_vec, ind) {
  vals <- unlist(lapply(list_vec, function(x) x[ind]))
  return(vals)
}

.wrap <- function(vec1, vec2) {
  out <- list()
  for(i in 1:length(vec1)) {
    out <- c(out, list(vec1[i], vec2[i]))
  }
  return(out)
}

pseudoCountNorm <- function(vals) {
  min.nz <- min(vals[vals!=0])
  vals[vals==0] <- min.nz
  return(vals)
}

minNAVal <- function(vals) {
  min.nz <- min(vals, na.rm=TRUE)
  vals[is.na(vals)] <- min.nz
  return(vals)
}
