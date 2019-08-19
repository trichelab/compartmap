#' Title
#'
#' @param obj 
#' @param genome 
#'
#' @return
#' @export
#'
#' @examples
filterOpenSea <- function(obj, genome = c("hg19", "hg38", "mm10", "mm9")) {
  #get the desired open sea loci given the genome
  genome <- match.arg(genome)
  openseas.genome <- switch(genome,
                            hg19=data("openSeas.hg19", package="compartmap"),
                            hg38=data("openSeas.hg38", package="compartmap"),
                            mm10=data("openSeas.mm10", package="compartmap"),
                            mm9=data("openSeas.mm9", package="compartmap"))
  #Subset by overlaps
  message("Filtering to open sea CpG loci...")
  obj.openseas <- subsetByOverlaps(obj, openseas.genome)
  return(obj.openseas)
}

maskArrays <- function(obj, genome = c("hg19", "hg38"),
                       maf = 0.1) {
  #mask arrays for bad actor probes and SNPs
}

getPMD <- function(obj, genome = c("hg19", "hg38")) {
  #return PMDs derived from Zhou et al. 2018 Nature Genetics
  
}

#' Title
#'
#' @param gr 
#' @param cutoff 
#'
#' @return
#' @export
#'
#' @examples
extractOpenClosed <- function(gr, cutoff = 0){
  #check for input to be GRanges
  if (!is(gr, "GRanges")) stop("Input needs to be a GRanges.")
  ifelse(gr$pc < cutoff, "open", "closed")
}

checkAssayType <- function(obj) {
  #helper function to check the class of an object
  is(obj, "SummarizedExperiment")
}

getAssayNames <- function(se) {
  #helper function to check the assay slot names
  names(assays(se))
}

cleanAssayRows <- function(se, rowmax = 0.5,
                           assay = c("array", "atac", "bisulfite")) {
  switch(assay,
         array = se[rowMeans(is.na(assays(se)$Beta)) < rowmax,],
         atac = se[rowMeans(is.na(assays(se)$counts)) < rowmax,],
         bisulfite = se[rowMeans(is.na(assays(se)$counts)) < rowmax,])
}

cleanAssayCols <- function(se, colmax = 0.8,
                           assay = c("array", "atac", "bisulfite")) {
  switch(assay,
         array = se[,colMeans(is.na(assays(se)$Beta)) < colmax],
         atac = se[,colMeans(is.na(assays(se)$counts)) < colmax],
         bisulfite = se[,colMeans(is.na(assays(se)$counts)) < colmax])
}

#' Helper function: squeezed logit
#'
#' @param x       a vector of values between 0 and 1 inclusive
#' @param sqz     the amount by which to 'squeeze', default is .000001
#'
#' @return        a vector of values between -Inf and +Inf
#'
#' @import        gtools
#'
#' @export 
flogit <- function(x, sqz=0.000001) {
  x[ which(x < sqz) ] <- sqz 
  x[ which(x > (1 - sqz)) ] <- (1 - sqz)
  logit(x)
}

#' Helper function: expanded expit
#'
#' @param x       a vector of values between -Inf and +Inf
#' @param sqz     the amount by which to 'squeeze', default is .000001
#'
#' @return        a vector of values between 0 and 1 inclusive
#'
#' @import        gtools
#'
#' @export 
fexpit <- function(x, sqz) {
  (((((inv.logit(x) * 2) - 1) / (1 - sqz)) + 1) / 2)
}
