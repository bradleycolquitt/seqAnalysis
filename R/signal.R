source("~/src/seqAnalysis/R/util.R")
library(parallel)
library(foreach)
library(doMC)
library(itertools)
registerDoMC(cores=6)

# Input vector of densities (e.g. H3K27ac density over enhancer region)
# Outputs Nx5 matrix, where N is number of valid dips found, and columns are
#   1) start of down slope
#   2) end of down slope
#   3) start of up slope
#   4) end of up slope
#   5) minimum value at floor (defined as between end of down and start of up slopes) for
#       dip ranking
find_min <- function(vals, window=5, height_min=5, floor_max=2, dip_min=2, dip_max=30) {
  
  # loess smooth values
  vals <- predict(loess(vals~c(1:length(vals)), span=0.05))
  diffs <- diff(vals)
  diffs.b <- ifelse(diffs>0, 1, -1)
  diffs.b[diffs==0] <- 0
  diffs.rle <- rle(diffs.b)
  
  # find continuous downs and ups
  diffs.down <- which(diffs.rle$lengths >= window & diffs.rle$values < 0)
  diffs.up <- which(diffs.rle$lengths >= window & diffs.rle$values > 0)
  
  # extract positions (both start and stop of up and down slopes)
  # returns 2xN matrix, where N is the number of slopes found
  diffs.down.pos <- sapply(diffs.down, function(x) c(sum(diffs.rle$lengths[1:(x-1)]), 
                                                     sum(diffs.rle$lengths[1:x])))
  diffs.up.pos <- sapply(diffs.up, function(x) c(sum(diffs.rle$lengths[1:(x-1)]),
                                                     sum(diffs.rle$lengths[1:x])))
  if (ncol(diffs.down.pos)==0 | ncol(diffs.up.pos)==0) {
    return(matrix(0, nrow=0, ncol=1))
  }
  
  # determine distances (from the end of the down slope to the start of the up slope)
  diffs.dist <- t(apply(diffs.down.pos, 2, function(down) apply(diffs.up.pos, 2, function(up) {
    up[1] - down[2]
  })))
  
  if (prod(dim(diffs.dist))==0) {
    return(matrix(0, nrow=0, ncol=1))
  }  
  
  # select pair of slopes within range
  valid_down_slopes <- list()
  valid_up_slopes <- list()
  for (i in 1:nrow(diffs.dist)) {
    row <- diffs.dist[i,]
    row.pos <- row>=0
    if (sum(row.pos) == 0) next
    min_val <- min(row[row.pos])
    if (min_val >= dip_min & min_val <= dip_max) {
      valid_down_slopes <- c(valid_down_slopes, i)
      valid_up_slopes <- c(valid_up_slopes, which(row==min_val))
    } 
  }

  if (length(valid_down_slopes)==1 | length(valid_up_slopes)==1) {
    #print("here")
    return(matrix(0, nrow=0, ncol=1))
  }  
  valid_slopes <- t(rbind(diffs.down.pos[,unlist(valid_down_slopes)], diffs.up.pos[,unlist(valid_up_slopes)]))
  #return(valid_slopes)
  #print(dim(valid_slopes))
  #print(prod(dim(valid_slopes)))
  if (prod(dim(valid_slopes))==0) {
    #print("here")
    return(matrix(0, nrow=0, ncol=1))
  } 

  # select dips within height range
  # extract minimum of tops of slopes and the minimum of the floor
  vals_select <- t(apply(valid_slopes, 1, function(x) c(min(vals[x[1]], vals[x[4]]), 
                                                        min(vals[x[2]:x[3]]))))

  valid_heights <- list()
  for (i in 1:nrow(vals_select)) {
    height <- abs(diff(vals_select[i,]))
    if (height >= height_min & vals_select[i,][2] <= floor_max) valid_heights <- c(valid_heights, i)
  }
  
  valid_heights <- unlist(valid_heights)
  valid_dips <- matrix(valid_slopes[valid_heights,], ncol=4)
  vals_select <- vals_select[valid_heights, 2]
  valid_dips <- cbind(valid_dips, vals_select)
  
  # valid_dips is now Nx5 matrix
  
  #return(valid_dips)
  # remove dips with common start or stop
  valid_dips <- valid_dips[!(duplicated(valid_dips[,2]) | duplicated(valid_dips[,3])),]
  valid_dips <- matrix(valid_dips, ncol=5, dimnames=list(c(), c("down_start",
                                                                "down_end",
                                                                "up_start",
                                                                "up_end",
                                                                "min_score")))
  return(valid_dips)
}  

filter_blank_dips <- function(dip_list) {
  ind <- lapply(dip_list, function(x) prod(dim(x)))
  return(dip_list[ind>0])
  
}
dip_to_bed <- function(bed, dip_matrix) {
  out_start_stops <- matrix(0, nrow=nrow(dip_matrix), ncol=3)
  for (i in 1:nrow(dip_matrix)) {
    out_start_stops[i,1] <- bed[dip_matrix[i,1],2]
    out_start_stops[i,2] <- bed[dip_matrix[i,4],3]
    out_start_stops[i,3] <- dip_matrix[i,5]
  }
  out_bed <- data.frame(bed[1,1], out_start_stops[,1], out_start_stops[,2], paste(bed[1,4], 1:nrow(dip_matrix), sep="."), out_start_stops[,3], bed[1,6])
  return(out_bed)
}


find_dips <- function(signal_matrix, bed, window=5, height_min=5, floor_max=2, dip_min=5, dip_max=30 ) {
  print("Finding dips...")
  dips <- foreach(x=iter(signal_matrix, by="row")) %dopar% find_min(x[1,], 
                                                                    window=window, 
                                                                    height_min=height_min, 
                                                                    floor_max=floor_max, 
                                                                    dip_min=dip_min, 
                                                                    dip_max=dip_max)
  
  names(dips) <- rownames(signal_matrix)
  dips.n <- filter_blank_dips(dips)
  print("Extracting bed records...")
  dips.bed <- mclapply(1:length(dips.n), function(x) dip_to_bed(bed[bed[,4] %in% names(dips.n)[x],], dips.n[[x]]), mc.cores=4, mc.preschedule=TRUE)
  dips.bed <- do.call("rbind", dips.bed)
  colnames(dips.bed) <- c("chr", "start", "stop", "id", "min_score", "strand")
  return(dips.bed)


}