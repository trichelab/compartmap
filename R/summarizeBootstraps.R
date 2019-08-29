#' Summarize the bootstrap compartment estimates and compute Agresti-Coull confidence intervals
#'
#' @name summarizeBootstraps
#'
#' @param boot.list List of bootstraps as GRanges objects 
#' @param est.ab The original compartment calls
#' @param q Confidence interval to compute (0.95 for 95% CI)
#' @param assay Type of assay we are working with
#'
#' @return A GRanges object with bootstraps summarized
#' @import SummarizedExperiment
#' @export
#'
#' @examples
#' 
summarizeBootstraps <- function(boot.list, est.ab, q = 0.95,
                                assay = c("array", "atac", "bisulfite")) {
  #go through the estimated A/B compartments and compute proportions from the boot.list
  est.ab$gr$score <- est.ab$pc
  message("Summarizing bootstraps.")
  #initialize open and closed counts to enumerate
  est.ab$gr$boot.open <- 0
  est.ab$gr$boot.closed <- 0
  #filter out failed compartment estimates
  boot.list <- removeEmptyBoots(boot.list)
  #summarize
  lapply(boot.list, function(b) {
    #add the pc to the granges object
    b$gr$score <- b$pc
    #convert to binary result for proportions
    #logic for methylation
    #open = eigen < 0
    #closed = eigen > 0
    #the logic is flipped for ATAC
    #open = eigen > 0
    #closed = eigen < 0
    if (assay == "atac") {
      b$gr$open <- ifelse(b$gr$score > 0, 1, 0)
      b$gr$closed <- ifelse(b$gr$score < 0, 1, 0)
    } else {
      b$gr$open <- ifelse(b$gr$score < 0, 1, 0)
      b$gr$closed <- ifelse(b$gr$score > 0, 1, 0)
    }
    #overlap by common intervals
    ol <- findOverlaps(b$gr, est.ab$gr)
    mcols(est.ab$gr)$boot.open[subjectHits(ol)] <<- mcols(est.ab$gr)$boot.open[subjectHits(ol)] + mcols(b$gr)$open[queryHits(ol)]
    mcols(est.ab$gr)$boot.closed[subjectHits(ol)] <<- mcols(est.ab$gr)$boot.closed[subjectHits(ol)] + mcols(b$gr)$closed[queryHits(ol)]
  })
  
  message("Computing Agresti-Coull 95% confidence intervals.")
  est.ab$gr$conf.est <- 0
  est.ab$gr$conf.est.upperCI <- 0
  est.ab$gr$conf.est.lowerCI <- 0
  lapply(1:length(est.ab$gr), function(e) {
    #compute intervals
    if (est.ab$gr[e,]$score > 0) {
      #assumes closed for ATAC
      est <- agrestiCoullCI(est.ab$gr[e,]$boot.closed,
                            est.ab$gr[e,]$boot.open,
                            q = 0.95)
      #this really isn't a good idea...
      est.ab$gr[e,]$conf.est.lowerCI <<- est[1]
      est.ab$gr[e,]$conf.est <<- est[2]
      est.ab$gr[e,]$conf.est.upperCI <<- est[3]
    }
    if (est.ab$gr[e,]$score < 0) {
      #assumes open for ATAC
      est <- agrestiCoullCI(est.ab$gr[e,]$boot.open,
                            est.ab$gr[e,]$boot.closed,
                            q = 0.95)
      #this really isn't a good idea...
      est.ab$gr[e,]$conf.est.lowerCI <<- est[1]
      est.ab$gr[e,]$conf.est <<- est[2]
      est.ab$gr[e,]$conf.est.upperCI <<- est[3]
    }
  })
  return(est.ab$gr)
}
