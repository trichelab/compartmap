#' Transform/normalize compartment calls using TF-IDF
#' 
#' @name transformTFIDF
#'
#' @param obj n x p input matrix (n = samples/cells; p = compartments)
#' @param scale.factor Scaling factor for the term-frequency (TF)
#'
#' @return A TF-IDF transformed matrix of the same dimensions as the input
#' @import Matrix
#' @export
#'
#' @examples
#' 
transformTFIDF <- function(obj, scale.factor = 1e5) {
  #this filters using TF-IDF on a *matrix* object
  if (!is(obj, "matrix")) stop("Input needs to be a matrix.")
  #the following code was modeled after or taken from:
  #http://andrewjohnhill.com/images/posts/2019-5-6-dimensionality-reduction-for-scatac-data/analysis.html
  #https://divingintogeneticsandgenomics.rbind.io/post/clustering-scatacseq-data-the-tf-idf-way/
  #binarize the matrix
  #this assumes n x p matrix (e.g. a wide matrix)
  #check and transpose as needed
  #input matrix is tall
  message("Input is a tall matrix. Transposing to wide.")
  if (dim(obj)[1] > dim(obj)[2]) obj <- t(obj)
  #make sparse
  obj.binary <- Matrix(.binarizeMatrix(t(obj)), sparse = TRUE)
  #compute term-frequency
  tf <- t(t(obj.binary) / Matrix::colSums(obj.binary))
  #scale
  tf@x <- log1p(tf@x * scale.factor)
  #inverse-document frequency smooth
  idf <- log(1 + ncol(obj.binary) / Matrix::rowSums(obj.binary))
  #transform
  tfidf <- .tfidf(tf, idf)
  #cast back to a matrix since things like UMAP don't like sparse matrices
  tfidf <- as.matrix(tfidf)
  return(t(tfidf))
}


#helper function
.binarizeMatrix <- function(obj) {
  #set positive values to 1 and negative to 0
  #open chromatin in atac or RNA is 1
  #closed chromatin in bisulfite or arrays is 1
  #just associate 1 with signal
  obj[obj > 0] <- 1
  obj[obj < 0] <- 0
  return(obj)
}

#helper function for TF-IDF transform
#modeled after http://andrewjohnhill.com/images/posts/2019-5-6-dimensionality-reduction-for-scatac-data/analysis.html
.tfidf <- function(tf, idf) {
  tf = t(tf)
  tf@x <- tf@x * rep.int(idf, diff(tf@p))
  tf = t(tf)
  return(tf)
}