#' Calculate Pearson correlations of smoothed eigenvectors
#' 
#' This function is used to generate a list x to be passed to getABSignal
#'
#' @param x      A list object from getCorMatrix
#' @param k    Value of k for smoothing (default = 2)
#' @param iter    Number of iterations for moving average smoothing (default = 2)
#' @param squeeze    Whether squeezing was used (implies Fisher's Z transformation)
#' 
#' @return    A list x to pass to getABSignal
#' 
#' @import    SummarizedExperiment
#' 
#' @export 
#' 
#' @examples 
#' 
#' library(SummarizedExperiment)
#' library(Homo.sapiens)
#' library(BiocSingular)
#' 
#' #Generate random genomic intervals of 1-1000 bp on chr1-22
#' #Modified from https://www.biostars.org/p/225520/
#' random_genomic_int <- data.frame(chr = rep("chr14", 100))
#' random_genomic_int$start <- apply(random_genomic_int, 1, function(x) { round(runif(1, 0, seqlengths(Homo.sapiens)[x][[1]]), 0) })
#' random_genomic_int$end <- random_genomic_int$start + runif(1, 1, 1000)
#' random_genomic_int$strand <- "*"
#' 
#' #Generate random counts
#' counts <- rnbinom(1000, 1.2, 0.4)
#' 
#' #Build random counts for 10 samples
#' count.mat <- matrix(sample(counts, nrow(random_genomic_int) * 10, replace = FALSE), ncol = 10)
#' colnames(count.mat) <- paste0("sample_", seq(1:10))
#' 
#' #Bin counts
#' bin.counts <- getBinMatrix(count.mat, makeGRangesFromDataFrame(random_genomic_int), chr = "chr14", genome = "hg19")
#' 
#' #Calculate correlations
#' bin.cor.counts <- getCorMatrix(bin.counts)
#' 
#' #Get A/B signal
#' absignal <- getABSignal(bin.cor.counts)

getABSignal <- function(x, k = 1, iter = 2, squeeze = FALSE){
  message("Calculating eigenvectors...")
  pc <- getSVD(x$binmat.cor, sing.vec = "right")
  if (squeeze) pc <- ifisherZ(pc)
  message("Smoothing with a k of ", k, " for ", iter, " iterations...")
  pc <- meanSmoother(pc, k=k, iter=iter)
  message("Done smoothing...")
  gr <- x$gr
  gr$pc <- pc
  gr$compartments <- extractOpenClosed(pc)
  return(gr)
}
