# returns "gene-level" gtf
process_gtf = function(gtf) {
  colnames(gtf) = c("chr", "build", "element", "start", "end", "x", "strand", "flag", "id")
  start = 0
  end = 0
  gtf_s = split(gtf, gtf[,9])
  gtf_coord = ldply(lapply(gtf_s, function(x) c(min(x[,4]), max(x[,5]))))

  gtf_un = gtf[!duplicated(gtf[,9]),]
  gtf_un[,3] = "gene"
  gtf_un[,4] = gtf_coord[,2]
  gtf_un[,5] = gtf_coord[,3]
  return(gtf_un)
 
}