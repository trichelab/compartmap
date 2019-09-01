#' Get the global means of a matrix
#'
#' @name getGlobalMeans
#'
#' @param obj Input SummarizedExperiment object
#' @param targets Column names or indices to indicate which samples to shrink towards
#'
#' @return A vector of global or targeted means
#' @import SummarizedExperiment
#' @export
#'
#' @examples
#' data("meth_array_450k_chr14", package = "compartmap")
#' array.global.means <- getGlobalMeans(array.data.chr14, assay = "array")
#'

getGlobalMeans <- function(obj, targets = NULL,
                           assay = c("array", "atac", "bisulfite")) {
  #match the assay arg
  assay <- match.arg(assay)

  #check the names of the assays
  if (!any(getAssayNames(obj) %in% c("Beta", "counts"))) {
    stop("The assay slot should contain either 'Beta' for arrays or 'counts' for atac/bisulfite.")
  }
  
  #get the global means
  #check if shrinkage targets are being used
  if (!is.null(targets)){
    stargets <- getShrinkageTargets(obj, targets)
    message("Using ", paste(shQuote(targets), collapse = ", "), " as shrinkage targets...")
    globalMean <- switch(assay,
                         array = matrix(rowMeans(flogit(assays(stargets)$Beta), na.rm=TRUE), ncol=1),
                         atac = matrix(rowMeans(assays(stargets)$counts, na.rm=TRUE), ncol=1),
                         bisulfite = matrix(rowMeans(assays(stargets)$counts, na.rm=TRUE), ncol=1))
  }
  else {
    globalMean <- switch(assay,
                         array = matrix(rowMeans(flogit(assays(obj)$Beta), na.rm=TRUE), ncol=1),
                         atac = matrix(rowMeans(assays(obj)$counts, na.rm=TRUE), ncol=1),
                         bisulfite = matrix(rowMeans(assays(obj)$counts, na.rm=TRUE), ncol=1))
  }
  colnames(globalMean) <- "globalMean"
  return(globalMean) 
}

#helper function for resampling
.getResampledGlobalMeans <- function(obj, targets = NULL,
                                     assay = c("array", "atac", "bisulfite")) {
  #match the assay arg
  assay <- match.arg(assay)
  
  #check the names of the assays
  if (!any(getAssayNames(obj) %in% c("Beta", "counts"))) {
    stop("The assay slot should contain either 'Beta' for arrays or 'counts' for atac/bisulfite.")
  }
  
  #get the global means
  #check if shrinkage targets are being used
  if (!is.null(targets)){
    stargets <- getShrinkageTargets(obj, targets)
    if (ncol(stargets) < 5) {
      message("Not enough samples to bootstrap using the shrinkage targets.")
      stop("Need at least 5 samples for this to be useful.")
      }
    globalMean <- switch(assay,
                         array = matrix(rowMeans(flogit(assays(stargets)$Beta), na.rm=TRUE), ncol=1),
                         atac = matrix(rowMeans(assays(stargets)$counts, na.rm=TRUE), ncol=1),
                         bisulfite = matrix(rowMeans(assays(stargets)$counts, na.rm=TRUE), ncol=1))
  }
  else {
    globalMean <- switch(assay,
                         array = matrix(rowMeans(flogit(assays(obj)$Beta), na.rm=TRUE), ncol=1),
                         atac = matrix(rowMeans(assays(obj)$counts, na.rm=TRUE), ncol=1),
                         bisulfite = matrix(rowMeans(assays(obj)$counts, na.rm=TRUE), ncol=1))
  }
  colnames(globalMean) <- "globalMean"
  return(globalMean) 
}
