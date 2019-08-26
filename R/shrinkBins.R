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
#' @param targets the column/sample/cell names to shrink towards
#' @param resample 
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

shrinkBins <- function(x, prior.means = NULL, chr = NULL,
                       res = 1e6, targets = NULL, resample = FALSE,
                       assay = c("array", "atac", "bisulfite"), ...) {
  #match the assay args
  assay <- match.arg(assay)
  
  #double check the obj class is compatible
  if (!checkAssayType(x)) stop("Input needs to be a SummarizedExperiment")
  
  #get the prior means
  if (is.null(prior.means)) {
    prior.means <- getGlobalMeans(obj=x, targets=targets, assay=assay)
  }
  #if we are resampling 
  if (is.null(prior.means) & resample) {
    prior.means <- .getResampledGlobalMeans(obj=x, targets=targets, assay=assay)
    }

  #bin the input
  #atac counts are summed
  if (assay == "atac") {
    bin.mat <- getBinMatrix(x=as.matrix(cbind(assays(x)$counts, prior.means)),
                            genloc=rowRanges(x), chr=chr, res=res, FUN=sum, ...)
  } else {
    #array and bisulfite values are taken as the median value
    bin.mat <- getBinMatrix(x=as.matrix(cbind(assays(x)$counts, prior.means)),
                            genloc=rowRanges(x), chr=chr, res=res, FUN=median, ...)
  }
  
  #shrink the bins using a James-Stein Estimator
  x.shrink <- apply(bin.mat$x, 1, function(r) {
    r.samps <- r[!names(r) %in% "globalMean"]
    r.prior.m <- r["globalMean"]
    if (!is.null(targets)) {
      if (length(r.samps[targets]) == 1) {
        stop("Cannot perform targeted bin-level shrinkage with one target sample.")
      }
      switch(assay,
             atac = .shrinkATAC(x=r.samps, prior=r.prior.m, targets=targets),
             array = .shrinkArray(x=r.samps, prior=r.prior.m, targets=targets),
             bisulfite = .shrinkBS(x=r.samps, prior=r.prior.m, targets=targets))
    }})
  
  return(list(gr=bin.mat$gr, x=x.shrink, gmeans=bin.mat$x[,"globalMean"]))
}

#helper functions for computing shrunken means
#shrink bins in (sc)ATAC-seq data
.shrinkATAC <- function(x, prior = NULL, offset = 0.0001,
                        targets = NULL) {
  if (!is.null(targets)) {
    C <- sd(log(x[targets] + offset))
  } else {
    C <- sd(log(x + offset))
  }
  prior.m <- log(prior + offset)
  #convert back to counts
  return(round(exp(prior.m + C*(log(x + offset) - prior.m))))
}

#shrink bins in methylation arrays
.shrinkArrays <- function(x, prior = NULL, targets = NULL) {
  if (!is.null(targets)) {
    C <- sd(flogit(x[targets]))
  } else {
    C <- sd(flogit(x))
  }
  prior.m <- flogit(prior)
  #convert back to beta values
  return(fexpit(prior.m + C*(flogit(x) - prior.m)))
}

#shrink bisulfite sequencing smoothed M-values
.shrinkBS <- function(x, prior = NULL, targets = NULL) {
  #assumes that M-values exist already
  if (!is.null(targets)) {
    C <- sd(x[targets])
  } else {
    C <- sd(x)
  }
  prior.m <- prior
  return(prior.m + C*(x - prior.m))
}
