#' @title Estimate A/B compartments from ATAC-seq data
#'
#' @description 
#' \code{getATACABsignal} returns estimated A/B compartments from ATAC-seq data.
#'
#' @param obj Input SummarizedExperiment object
#' @param res Compartment resolution in bp
#' @param parallel Whether to run samples in parallel
#' @param chr What chromosome to work on (leave as NULL to run on all chromosomes)
#' @param targets Samples/cells to shrink towards
#' @param cores How many cores to use when running samples in parallel
#' @param bootstrap Whether we should perform bootstrapping of inferred compartments
#' @param num.bootstraps How many bootstraps to run
#' @param genome What genome to work on ("hg19", "hg38", "mm9", "mm10")
#' @param other Another arbitrary genome to compute compartments on
#' @param group Whether to treat this as a group set of samples
#'
#' @return A RaggedExperiment of inferred compartments
#' @import SummarizedExperiment
#' @import parallel
#' @import RaggedExperiment
#' @export
#' @examples
#' data("groupATAC_raw_filtered_chr14", package = "compartmap")
#' atac_compartments <- getATACABsignal(filtered.data.chr14, parallel=F, chr="chr14", bootstrap=F, genome="hg19")

scCompartments <- function(obj, res = 1e6, parallel = TRUE, chr = NULL,
                           targets = NULL, cores = 2,
                           bootstrap = TRUE, num.bootstraps = 1000,
                           genome = c("hg19", "hg38", "mm9", "mm10"),
                           group = FALSE, assay = c("atac", "bisulfite", "rna")) {
  #this is a _complete_ rebuild of getATACABsignal()
  #the intended purpose is to refactor and streamline the codebase to make this more usable
  
  #gather the chromosomes we are working on
  if (is.null(chr)) {
    message("Assuming we want to process all chromosomes.")
    #get what chromosomes we want
    chr <- getChrs(obj)
  }
  
  #get the column names
  if (is.null(colnames(obj))) stop("colnames needs to be sample names.")
  columns <- colnames(obj)
  names(columns) <- columns
  
  #precompute global means
  prior.means <- getGlobalMeans(obj = obj, targets = targets, assay = "atac")
  
  if (bootstrap) {
    message("Pre-computing the bootstrap global means.")
    bmeans <- precomputeBootstrapMeans(obj = obj, targets = targets, num.bootstraps = num.bootstraps,
                                       assay = "atac", parallel = parallel, num.cores = cores)
  }
  
  #this is where we should shrink bins a priori because it's too slow otherwise
  prior.means.shrink.bins <- lapply(chr, function(c) {
    message("Shrinking bins for chromosome ", c)
    pmeans <- as(prior.means, "GRanges")
    pmeans <- keepSeqlevels(pmeans, chr, pruning.mode = "coarse")
    #go back to a matrix
    prior.means <- as(pmeans, "matrix")
    colnames(prior.means) <- "globalMean"
    #get the shrunken bins
    obj.sub <- obj[,colnames(obj)]
    obj.bins <- shrinkBins(obj.sub, obj, prior.means = pmeans, chr = c,
                           res = res, targets = targets, assay = "atac",
                           genome = genome)
    return(obj.bins)
  })
  
  #set the names to chromosomes
  #each element in the list should correspond the chromosomes analyzed
  names(prior.means.shrink.bins) <- chr
  
  #need to do this for bootstraps too... this is trickier
  #nested lists? yuck.
  if (bootstrap) {
    bootstrap.shrink.bins <- lapply(1:num.bootstraps, function(b) {
      prior.means.shrink.bins.boot <- lapply(chr, function(c) {
        message("Shrinking bootstrap bins for chromosome ", c)
        bmeans <- as(bmeans, "GRanges")
        bmeans <- keepSeqlevels(bmeans, c, pruning.mode = "coarse")
        #go back to a matrix
        bmeans <- as(bmeans, "matrix")
        b.sub <- as.matrix(bmeans[,b])
        colnames(b.sub) <- "globalMean"
        #get the shrunken bins
        #this is a kludgy catch...
        obj.sub <- obj[,colnames(obj)]
        obj.bins <- shrinkBins(obj.sub, obj, prior.means = b.sub, chr = c,
                               res = res, targets = targets, assay = "atac",
                               genome = genome)
        return(obj.bins)
      })
      #we know have a bootstrapped shrunken set of bins!
      names(prior.means.shrink.bins.boot) <- paste0("bootstrap_", b) 
    })
  }
  
  #pass on the prior.means as the pre-shrunken bins to the function
  #pass on the bootstrap.means as the list of pre-shrunken bins to the function
  return()
}