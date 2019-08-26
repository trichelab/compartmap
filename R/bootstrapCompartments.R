bootstrapCompartments <- function(se.obj, bootstrap.samples = 1000, with.ci = FALSE) {
  #function for nonparametric bootstrap of compartments to compute 95% CI
  
  #check input
  if (!is(se.obj, "RangedSummarizedExperiment") | !is(se.obj, "SummarizedExperiment")) {
    stop("Input needs to be a (Ranged)SummarizedExperiment object.")
  }
  
  #get the proportion of samples to take everytime
  cells.to.sample <- ceiling(ncol(se.obj) * sample.prop)
  if (cells.to.sample == ncol(se.obj)) stop("Somehow we are trying to sample all columns of the matrix. Decrease the sample.prop.")
  if (cells.to.sample == 1) stop("Too few samples to use for bootstrapping. Increase the sample.prop.")
  
  #resample the global means with replacement
  
  return()
}

#A helper function to reconstruct the GRanges object for binning if not supplied from the input matrix
#These *should* be the rownames of the input matrix
.buildRanges <- function(obj) {
  if (!is.null(rownames(obj))) {
    message("Constructing rowRanges from the input matrix...")
    snames <- sapply(strsplit(rownames(obj), ":"), `[`, 1)
    loc <- as.numeric(sapply(strsplit(rownames(obj), ":"), `[`, 2))
    myranges <- GRanges(seqnames = snames,
                        strand = "*",
                        ranges = IRanges(start = loc, width = 1))
  }
  else (stop("Could not reconstruct the rowRanges from the input matrix..."))
  return(myranges)
}

#helper function to re-sample
#this was inspired by https://github.com/sgibb/bootstrap/blob/master/R/helper-functions.R
.resampleMatrix <- function(x, size=ncol(x)) {
  samp.to.select <- sample.int(ncol(x), size=size, replace=TRUE)
  return(x[, samp.to.select])
}