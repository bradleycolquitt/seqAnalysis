## Intersect/setdiff feature regions of two different samples
MP.feature.compare <- function(sample1, sample2, method="intersect") {
  s1.name <- sample1$name
  s2.name <- sample2$name
  if (method == "intersect") result <- s1.name[s1.name%in%s2.name]
  return(result)
}


MP.tabulateSig <- function() {
  samples <- c(samples.hmc, samples.mc)
  out.1 <- foreach(sample=samples, .verbose=FALSE,
                   .inorder=TRUE, .combine="rbind") %dopar% {
    sample.path <- paste(mp.path, sample, sep="/")
    features <- list.files(sample.path, pattern="thresh")
    out.2 <- foreach(feature=features, .combine="rbind") %do% {
      feature.path <- paste(sample.path, feature, sep="/")
      sigs <- list.files(feature.path)
      if (length(sigs) > 0) {
        out.3 <- foreach(i=1:length(sigs), .combine="rbind") %do% {
          name <- unlist(strsplit(sigs[i], "_"))
          out.tmp <- c(name[1], name[2], name[3])
          lines <- length(count.fields(paste(feature.path, sigs[i], sep="/")))
          out.tmp <- c(out.tmp, lines)
          return(out.tmp)
        }
        return(cbind(rep(feature, times=nrow(out.3)), out.3))
      }
    }
    return(cbind(rep(sample, times=nrow(out.2)), out.2))
  }
  colnames(out.1) <- c("samples", "feature",  "FDR", "ratio", "direction", "number")
  out.1 <- out.1[-grep(".bed", out.1[,'direction']),]
  out.1 <- out.1[-grep("test", out.1[,'FDR']),]
  out.1 <- transform(out.1, samples <- factor(samples, levels=unique(samples)),
                     feature <- factor(feature, levels=unique(feature)),
                     FDR <- factor(FDR, levels=unique(FDR)),
                     ratio <- factor(ratio, levels=unique(ratio)),
                     direction <- factor(direction, levels=unique(direction)),
                     number <- as.numeric(number))
  write.table(out.1, file=paste(mp.path, "feature_summaries/diff_enrich_tab.txt", sep="/"),
              quote=FALSE, sep="\t", row.names=FALSE)
  return(out.1)
}

MP.intersectFeature <- function(mp, feature) {
  mp.data <- read.delim(mp)
  feature.data <- read.delim(feature)
  mp.split <- split(mp.split, mp.split$chr)
  feature.split <- read.delim(feature, feature[,1])
  matches <- foreach(mp=mp.split, feature=feature.split) %dopar% {
    foreach (feat=isplitRows(feature, chunk=1)) %do% {
      return(mp[mp$start >= feat[,2] & mp$stop <= feat[,3]])
    }
  }
  return(matches)
                         
}


## Extract feature regions by thresholded key values
MP.thresh <- function(sample, feature, key="ams_A", thresh, above=TRUE) {
  data <- getFeatureFile(sample, feature)
  data.thresh <- dataThresh(data, key, thresh, above)
  return(data.thresh)
}

MP.thresh.byQ <- function(vals, N, above=TRUE) {
  q <- quantile(vals, probs=(length(vals)-N) / length(vals))
  out <- vals >= q
  return(out)
}


MP.plotMeans <- function(vals, select=NULL, columns=c(1,6), fname=NULL) {
  if (is.null(fname)) {X11("", 10, 10)
  } else {png(file=paste(cluster.path, "plots", fname, sep="/"), width=720, height=720)}  
  if (!is.null(select)) vals <- vals[match(select, rownames(vals)),]
  plot(1, 1, type="n", xlim=c(0, nrow(vals) + 1), ylim=c(200, 700))
  for(i in seq(columns[1], columns[2])) {
     points(vals[,i], pch=20, col=mean.cols[i])
  }
}


## Calculate mean AMS levels for given sample and feature
MP.computeMean <- function(sample, feature) {
  data <- getFeatureFile(sample, feature)
  val <- mean(data$ams_A)
  return(val)
}

MP.computeSD <- function(sample, feature) {
  data <- getFeatureFile(sample, feature)
  val <- sd(data$ams_A)
  return(val)
}

MP.computeStat.all <- function(margin=2, FUN=mean, test=FALSE, fname=NULL) {
  files <- list.files(feature.path)
  vals <- foreach(file=files, .combine="rbind", .inorder=TRUE) %do% {
    print(file)
    cells.hmc <- getFeatureFile(feature=file, sample=samples.hmc)
    cells.mc <-  getFeatureFile(feature=file, sample=samples.mc)
    cells.vals <- removeInf(log(cbind(cells.hmc[c('ams_A', 'ams_B')], cells.mc[,c('ams_A', 'ams_B')]),2))
    if (test) {
      comp <- list(c(1,2), c(3,4))
      cells.compute <- lapply(comp, function(x) {
                              t.test(cells.vals[,x[1]], cells.vals[,x[2]])$p.value
                            })
    } else {cells.compute <- apply(cells.vals, margin, FUN)
    }
    if (margin == 1) cells.mean <- mean(cells.mean)
    return(cells.compute)
  }
  rownames(vals) <- files
  if (!is.null(fname)) write.table(vals, file=paste(mp.path, "feature_summaries", fname, sep="/"), quote=FALSE, sep="\t")
  return(vals)
}

.cv <- function(x) {
  return(sd(x)/mean(x))
}
