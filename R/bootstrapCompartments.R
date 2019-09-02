#' Non-parametric bootstrapping of compartments and summarization of bootstraps/compute confidence intervals
#'
#' @name bootstrapCompartments
#'
#' @param obj List object of computed compartments for a sample with 'pc' and 'gr' as elements
#' @param original.obj The original, full input SummarizedExperiment of all samples/cells
#' @param bootstrap.samples How many bootstraps to run
#' @param chr Which chromosome to operate on
#' @param assay What sort of assay are we working on
#' @param parallel Whether to run the bootstrapping in parallel
#' @param cores How many cores to use for parallel processing
#' @param targets Targets to shrink towards
#' @param res The compartment resolution
#' @param genome What genome are we working on
#' @param q What sort of confidence intervals are we computing (e.g. 0.95 for 95\% CI)
#' @param svd The original compartment calls as a GRanges object
#' @param group Whether this is for group-level inference
#'
#' @return Compartment estimates with summarized bootstraps and confidence intervals
#' @import parallel
#' @import SummarizedExperiment
#' @export
#'
#' @examples
#' 
bootstrapCompartments <- function(obj, original.obj, bootstrap.samples = 1000,
                                  chr = "chr14", assay = c("array", "atac", "bisulfite"),
                                  parallel = TRUE, cores = 2, targets = NULL, res = 1e6,
                                  genome = c("hg19", "hg38", "mm9", "mm10"), q = 0.95,
                                  svd = NULL, group = FALSE) {
  #function for nonparametric bootstrap of compartments and compute 95% CIs
  #check input
  #match the assay args
  assay <- match.arg(assay)

  #double check the obj class is compatible
  if (!checkAssayType(original.obj)) stop("Input needs to be a SummarizedExperiment")
  
  #check the names of the assays
  if (!any(getAssayNames(original.obj) %in% c("Beta", "counts"))) {
    stop("The assay slot should contain either 'Beta' for arrays or 'counts' for atac/bisulfite.")
  }
  
  #if we are using targeted means
  if (!is.null(targets)) original.obj <- original.obj[,targets]
  if (ncol(original.obj) < 6) stop("We need more than 5 samples to bootstrap with for the results to be meaningful.")
  
  if (!parallel) {
    message("Not bootstrapping in parallel will take a long time...")
    #bootstrap and recompute compartments
    resamp.compartments <- lapply(1:bootstrap.samples, function(b) {
      #resample the global means with replacement
      message("Working on bootstrap ", b)
      resamp.mat <- switch(assay,
                           array = .resampleMatrix(assays(original.obj)$Beta),
                           atac = .resampleMatrix(assays(original.obj)$counts),
                           bisulfite = .resampleMatrix(assays(original.obj)$counts))
      #turn back into SummarizedExperiment
      resamp.se <- switch(assay,
                          array = SummarizedExperiment(assays=SimpleList(Beta=resamp.mat),
                                                       rowRanges = rowRanges(original.obj)),
                          atac = SummarizedExperiment(assays=SimpleList(counts=resamp.mat),
                                                      rowRanges = rowRanges(original.obj)),
                          bisulfite = SummarizedExperiment(assays=SimpleList(counts=resamp.map),
                                                           rowRanges = rowRanges(original.obj)))
      #get the shrunken bins with new global mean
      s.bins <- shrinkBins(obj, original.obj, prior.means = getGlobalMeans(resamp.se, targets = targets, assay = assay),
                           chr = chr, res = res, assay = assay, genome = genome)
      if (group) cor.bins <- getCorMatrix(s.bins, squeeze = FALSE)
      if (isFALSE(group)) cor.bins <- getCorMatrix(s.bins, squeeze = TRUE)
      #Stupid check for perfect correlation with global mean
      if (any(is.na(cor.bins$binmat.cor))) {
        absig <- matrix(rep(NA, nrow(cor.bins$binmat.cor)))
      } else {
        absig <- getABSignal(cor.bins, assay = assay)
      }
      return(absig)
    })
  } else {
    message("Bootstrapping in parallel with ", cores, " cores.")
    #bootstrap and recompute compartments
    resamp.compartments <- mclapply(1:bootstrap.samples, function(b) {
      #resample the global means with replacement
      resamp.mat <- switch(assay,
                           array = .resampleMatrix(assays(original.obj)$Beta),
                           atac = .resampleMatrix(assays(original.obj)$counts),
                           bisulfite = .resampleMatrix(assays(original.obj)$counts))
      #turn back into SummarizedExperiment
      resamp.se <- SummarizedExperiment(assays=SimpleList(Beta=resamp.mat),
                                        rowRanges = rowRanges(original.obj))
      #get the shrunken bins with new global mean
      s.bins <- shrinkBins(obj, original.obj, prior.means = getGlobalMeans(resamp.se, targets = targets, assay = assay),
                           chr = chr, res = res, assay = assay, genome = genome)
      cor.bins <- getCorMatrix(s.bins, squeeze = TRUE)
      #Stupid check for perfect correlation with global mean
      if (any(is.na(cor.bins$binmat.cor))) {
        absig <- matrix(rep(NA, nrow(cor.bins$binmat.cor)))
      } else {
        absig <- getABSignal(cor.bins, assay = assay)
      }
      return(absig)
    }, mc.cores = cores)
  }
  
  #summarize the bootstraps and compute confidence intervals
  resamp.compartments <- summarizeBootstraps(resamp.compartments, svd,
                                             q = q, assay = assay)
  return(resamp.compartments)
}

#helper function to re-sample
#this was inspired by https://github.com/sgibb/bootstrap/blob/master/R/helper-functions.R
.resampleMatrix <- function(x, size=ncol(x)) {
  samp.to.select <- sample.int(ncol(x), size=size, replace=TRUE)
  return(x[, samp.to.select])
}