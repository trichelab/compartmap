#' @title Estimate A/B compartments
#'
#' @description 
#' \code{getCompartments} returns estimated A/B compartments from ATAC-seq and methylation array data
#'
#' @details 
#' This is a wrapper function to perform A/B compartment inference. Compartmentalizer implements a Stein estimator to shrink per-sample compartment estimates towards a global mean. The expected input for this function can be generated using packages like SeSAMe and ATACseeker.
#'
#' @param obj The object with which to perform compartment inference
#' @param type The type of data that obj represents (e.g. atac or array)
#' @param res Resolution of compartments in base pairs (default is 1e6)
#' @param parallel Should the estimates be done in parallel (default is FALSE)
#' @param chrs Chromosomes to operate on (can be individual chromosomes, a list of chromosomes, or all)
#' @param genome Genome to use (default is hg19)
#' @param targets Specify samples to use as shrinkage targets
#' @param cores Specify the number of cores to use when running in parallel
#' @param run_examples Whether to run ATAC-seq and 450k example analysis
#' @param ... Other parameters to pass to internal functions
#'
#' @return A p x n matrix (samples as columns and compartments as rows) to pass to embed_compartments
#' 
#' @import gtools 
#' @import parallel
#' @import Homo.sapiens
#' @import minfi
#' @import GenomicRanges
#' 
#' @export
#'
#' @examples
#' 
#' library(GenomicRanges)
#' library(SummarizedExperiment)
#' library(Homo.sapiens)
#' 
#' #ATAC-seq data
#' data(bulkATAC_raw_filtered_chr14, package = "compartmap")
#' atac_compartments <- getCompartments(filtered.data.chr14, type = "atac", parallel = FALSE, chrs = "chr14")
#' \dontrun{
#' #450k data
#' data(meth_array_450k_chr14, package = "compartmap")
#' array_compartments <- getCompartments(array.data.chr14, type = "array", parallel = FALSE, chrs = "chr14")}

getCompartments <- function(obj, type = c("atac", "array"), res = 1e6, parallel = FALSE,
                             chrs = "chr1", genome = "hg19", targets = NULL,
                             cores = 1, run_examples = FALSE, ...) {

  # short circuit type checking if just testing
  if (run_examples) return(.run_examples())

  # Perform initial check the input data type
  type <- match.arg(type)
  typeclass <- c(atac="RangedSummarizedExperiment", array="GenomicRatioSet")
  typepackage <- c(atac="ATACseeker", array="minfi (and/or sesame)")
  if (!is(obj, typeclass[type])) {
    message(type, " data must be supplied as a ", typeclass[type], " object.")
    stop(paste("You can generate one with the", typepackage[type], "package."))
  }
  
  #Pre-check the chromosomes to be analyzed
  allchrs <- (length(chrs) == 1 & chrs == "all") 
  if (allchrs) {
    chrs <- seqlevels(obj)
    message("Mapping all chromosomes...")
  } else {
    message("Mapping chromosome", ifelse(length(chrs) > 1, "s ", " "), 
            paste(shQuote(chrs), collapse=", "))
  }

  #set the number of cores if running in parallel
  if (parallel) options(mc.cores = cores)
  if (parallel & cores == 1) options(mc.cores = detectCores()/2)
  
  # call the appropriate function
  switch(type,
         atac=getATACABsignal(obj=obj, res=res, parallel=parallel, 
                              allchrs=allchrs, chr=chrs, targets=targets, ...),
         array=getArrayABsignal(obj=obj, res=res, parallel=parallel, 
                                allchrs=allchrs, chr=chrs, targets=targets,...))
  
}


#' Example ATAC-seq data for compartmap
#' 
#' This data was generated using the data from the reference via bwa mem
#' and pre-processing the data using the ATACseeker package.
#' 
#' @name filtered.data.chr14
#' @docType data
#' @author Benjamin K Johnson \email{ben.johnson@vai.org}
#' @references \url{https://trace.ncbi.nlm.nih.gov/Traces/sra/?study=SRP082417}
#' @keywords data
#' @usage data(bulkATAC_raw_filtered_chr14, package = "compartmap")
NULL

#' Example Illumina 450k methylation array data for compartmap
#' 
#' This data was generated using the data from the reference via the
#' sesamize function from the SeSAMe package.
#' 
#' @name array.data.chr14
#' @docType data
#' @author Benjamin K Johnson \email{ben.johnson@vai.org}
#' @references \url{https://f1000research.com/articles/5-1281/v3}
#' @keywords data
#' @usage data(meth_array_450k_chr14, package = "compartmap")
NULL

# helper fn
.run_examples <- function() { 
    message("Running ATAC-seq example data...")
    data(bulkATAC_raw_filtered_chr14, package = "compartmap")
    atac_compartments <- getCompartments(filtered.data.chr14, type = "atac", parallel = FALSE, chrs = "chr14")
    message("Done!")
    message("Running 450k example data...")
    data(meth_array_450k_chr14, package = "compartmap")
    array_compartments <- getCompartments(array.data.chr14, type = "array", parallel = FALSE, chrs = "chr14")
    message("Done!")
    return(list(atac = atac_compartments, array = array_compartments))
}
