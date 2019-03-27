#' @title Employ an eBayes shrinkage approach for bin-level estimates for A/B inference
#'
#' @description 
#' \code{shrinkBins} returns shrunken bin-level estimates
#'
#' @details This function computes shrunken bin-level estimates using a James-Stein estimator, reformulated as an eBayes procedure
#' 
#' @param x a RangedSummarizedExperiment
#' @param prior.means the means of the bin-level prior distribution
#' @param assay which assay this is
#' @param chr the chromosome to operate on
#' @param res resolution to perform the binning
#' @param ... other parameters to pass to the internal functions
#'
#' @return A list object to pass to getCorMatrix
#' 
#' @import GenomicRanges
#' @import SummarizedExperiment
#' 
#' @export
#'
#' @examples
#' 

shrinkBins <- function(x, prior.means = NULL, assay = c("atac", "array"),
                       chr = NULL, res = 1e6, ...) {
  
  #determine the input assay type
  switch(assay,
         atac = shrinkBinsATAC(x=x, prior.means=prior.means, chr=chr, res=res, ...),
         array = shrinkBinsArray(x=x, prior.means=prior.means, chr=chr, res=res, ...))
  
}

#' @title Employ an eBayes shrinkage approach for bin-level estimates for A/B inference in ATAC-seq
#'
#' @description 
#' \code{shrinkBinsATAC} returns shrunken bin-level estimates from ATAC-seq
#'
#' @details This function computes shrunken bin-level estimates using a James-Stein estimator, reformulated as an eBayes procedure
#' 
#' @param x a RangedSummarizedExperiment
#' @param prior.means the means of the bin-level prior distribution
#' @param chr the chromosome to operate on
#' @param res resolution to perform the binning
#' @param targets the subset of samples to use for targeted shrinkage
#' @param ... other parameters to pass to the internal functions
#'
#' @return A list object to pass to getCorMatrix
#' 
#' @import GenomicRanges
#' @import SummarizedExperiment
#' 
#' @export
#'
#' @examples
#' 

shrinkBinsATAC <- function(x, prior.means = NULL, chr = NULL, res = 1e6, targets = NULL, ...) {
  #get the prior means
  if (is.null(prior.means) & is.null(targets)) prior.means <- .getGlobalMeansATAC(x)
  if (is.null(prior.means) & !is.null(targets)) prior.means <- .getGlobalMeansATAC(x, targets)

  #bin the input
  if (!is(x, "RangedSummarizedExperiment")) stop("The input doesn't look like a RangedSummarizedExperiment")
  bin.mat <- getBinMatrix(x=as.matrix(cbind(assays(x)$counts, prior.means)), genloc=rowRanges(x), chr=chr, res=res, FUN=sum, ...)
  
  #shrink the bins using a James-Stein Estimator
  if (is.null(targets)) {
    x.shrink <- apply(bin.mat$x, 1, function(r) {
      r.samps <- r[!names(r) %in% "globalMean"]
      r.prior.m <- r["globalMean"]
      C <- sd(log(r.samps + 0.0001))
      r.prior.m <- log(r.prior.m + 0.0001)
      return(round(exp(r.prior.m + C*(log(r.samps + 0.0001) - r.prior.m))))
      })
    } else {
      #assumes targeted shrinkage
      x.shrink <- apply(bin.mat$x, 1, function(r) {
        r.samps <- r[!names(r) %in% "globalMean"]
        r.prior.m <- r["globalMean"]
        #if the prior mean is targeted, we need to account for that with the shrinkage parameter C
        if (length(r.samps[targets]) == 1) {
          message("Cannot perform targeted bin-level shrinkage with one target sample.")
          stop("Re-run without targets or bin.shrink = FALSE")
          }
        C <- sd(log(r.samps[targets] + 0.0001))
        r.prior.m <- log(r.prior.m + 0.0001)
        return(round(exp(r.prior.m + C*(log(r.samps + 0.0001) - r.prior.m))))
      })
    }
  
  return(list(gr=bin.mat$gr, x=x.shrink, gmeans=bin.mat$x[,"globalMean"]))
}

#' @title Employ an eBayes shrinkage approach for bin-level estimates for A/B inference in methylation arrays
#'
#' @description 
#' \code{shrinkBinsATAC} returns shrunken bin-level estimates from methylation arrays
#'
#' @details This function computes shrunken bin-level estimates using a James-Stein estimator, reformulated as an eBayes procedure
#' 
#' @param x a RangedSummarizedExperiment
#' @param prior.means the means of the bin-level prior distribution
#' @param chr the chromosome to operate on
#' @param res resolution to perform the binning
#' @param ... other parameters to pass to the internal functions
#'
#' @return A list object to pass to getCorMatrix
#' 
#' @import GenomicRanges
#' @import SummarizedExperiment
#' 
#' @export
#'
#' @examples
#' 

shrinkBinsArray <- function(x, prior.means = NULL, chr = NULL, res = 1e6, ...) {
  return(stop("Array bin-level shrinkage is not implemented yet..."))
}
