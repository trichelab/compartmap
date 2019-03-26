#' @title Estimate A/B compartments from ATAC-seq data
#'
#' @description 
#' \code{getATACABsignal} returns estimated A/B compartments from methylation array data.
#'
#' @details This function estimates A/B compartments shrinking towards a global mean of targets or across samples
#' 
#'
#' @param obj Input GenomicRatioSet object 
#' @param res Compartment resolution (in bp)
#' @param parallel Should the inference be done in parallel?
#' @param allchrs Whether all autosomes should be used for A/B inference
#' @param chr Specify a chromosome to analyze
#' @param targets Specify samples as shrinkage targets
#' @param ... Additional arguments
#'
#' @return A p x n matrix (samples as columns and compartments as rows) of compartments
#' @import GenomicRanges
#' @import SummarizedExperiment
#' @import Homo.sapiens
#' @export
#' 
#' @examples 
#' library(GenomicRanges)
#' library(SummarizedExperiment)
#' library(Homo.sapiens)
#' 
#' data(bulkATAC_raw_filtered_chr14, package = "compartmap")
#' atac_compartments <- getATACABsignal(filtered.data.chr14, chr = "chr14", genome = "hg19")

getATACABsignal <- function(obj, res=1e6, parallel=FALSE, allchrs=FALSE, chr = NULL, targets = NULL, ...) {
  globalMeanSet <- .getGlobalMeansATAC(obj, targets)
  columns <- colnames(obj)
  names(columns) <- columns 
  
  getComp <- .getPaired
  if (allchrs == TRUE) getComp <- .getPairedAllChrs
  
  if (parallel) {
    #this should be set earlier...
    #options(mc.cores=detectCores()/2) # RAM blows up otherwise 
    do.call(cbind, 
            mclapply(columns,getComp,obj=obj,globalMeanSet=globalMeanSet,chr=chr,targets=targets,res=res,...))
  } else { 
    do.call(cbind, 
            lapply(columns,getComp,obj=obj,globalMeanSet=globalMeanSet,chr=chr,targets=targets,res=res,...))
  } 
}

.getGlobalMeansATAC <- function(obj, targets = NULL) {
  if (!is.null(targets)){
    stargets <- .getShrinkageTargets(obj, targets)
    message("Using ", paste(shQuote(targets), collapse = ", "), " as shrinkage targets...")
    meanBeta <- matrix(rowMeans(assay(stargets), na.rm=TRUE), ncol=1)
  }
  else (meanBeta <- matrix(rowMeans(assay(obj), na.rm=TRUE), ncol=1))
  colnames(meanBeta) <- "globalMean"
  return(meanBeta) 
}

.getShrinkageTargets <- function(obj, group) {
  if (all(group %in% colnames(obj))) stargets.obj <- obj[,group]
  else (stop("Could not find ", group, " in the colnames of the input matrix..."))
  return(stargets.obj)
}

.getPaired <- function(column, obj, globalMeanSet=NULL, res=1e6, chr=NULL, ...) {
  message("Computing shrunken compartment eigenscores for ", column, "...") 
  if(is.null(globalMeanSet)) globalMeanSet <- .getGlobalMeansATAC(obj)
  binmat <- getBinMatrix(as.matrix(cbind(assay(obj)[,column], globalMeanSet)), rowRanges(obj), chr = chr, res = res)
  cormat <- getCorMatrix(binmat, squeeze = TRUE)
  #Stupid check for perfect correlation with global mean
  if (any(is.na(cormat$binmat.cor))) {
    absig <- matrix(rep(NA, nrow(cormat$binmat.cor)))
  }
  else {
    absig <- getABSignal(cormat, squeeze = FALSE)
    absig.pc <- absig$pc
    names(absig.pc) <- as.character(granges(absig))
    absig <- absig.pc
  }
  return(absig)
}

.getPairedAllChrs <- function(column, obj, globalMeanSet=NULL, res=1e6, chr = NULL, ...) {
  if (is.null(globalMeanSet)) globalMeanSet <- .getGlobalMeansATAC(obj)
  #this isn't what we want to do...
  chrs <- paste0("chr", c(seq(1,22)))
  names(chrs) <- chrs
  getPairedChr <- function(chr) { 
    message("Computing shrunken eigenscores for ", column, " on ", chr, "...") 
    binmat <- getBinMatrix(cbind(assay(obj)[,column], globalMeanSet), rowRanges(obj), chr = chr, res = res)
    cormat <- getCorMatrix(binmat, squeeze = TRUE)
    #Stupid check for perfect correlation with global mean
    if (any(is.na(cormat$binmat.cor))) {
      absig <- matrix(rep(NA, nrow(cormat$binmat.cor)))
    }
    else {
      absig <- getABSignal(cormat, squeeze = FALSE)
      absig.pc <- absig$pc
      names(absig.pc) <- as.character(granges(absig))
      absig <- absig.pc
    }
    
    return(absig)
  }
  unlist(lapply(chrs, getPairedChr))
}
