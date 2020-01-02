#' Filter and convert/cluster compartment calls for embedding and dimension reduction
#'
#' @name embeddingPrep
#'
#' @param obj Input either a RaggedExperiment or output from fixCompartments
#' @param filter Whether to filter compartments
#' @param min.conf Minimum confidence estimate to use when filtering
#' @param min.eigen Minimum absolute eigenvalue to use when filtering
#' @param tf.idf Whether to TF-IDF transform the data
#' @param cluster Whether to identify clusters in the data
#' @param k How many neighbors when identifying clusters
#' @param dims How many dimensions to use for identifying clusters
#'
#' @return A prepared matrix of list of compartment calls for embedding and dimension reduction
#' @import RaggedExperiment
#' @import SummarizedExperiment
#' @import Matrix
#' @import BiocSingular
#' @import RANN
#' @import igraph
#' @export
#'
#' @examples
#' 
embeddingPrep <- function(obj, filter = TRUE,
                          min.conf = 0.7, min.eigen = 0.02,
                          tf.idf = FALSE, cluster = TRUE,
                          k = 10, dims = NULL) {
  
  #check and see if we have just one sample and stop if we do
  if (is(obj, "GRanges")) stop("Need at least 2 samples.")
  
  #fix up the compartments
  if (is(obj, "RaggedExperiment")) obj <- fixCompartments(obj, min.conf = min.conf)

  #filter compartments down
  if (isTRUE(filter)) {
    obj <- filterCompartments(obj,
                              min.conf = min.conf,
                              min.eigen = min.eigen)
  }
  
  #set a seed
  set.seed(1000)
  #turn into a matrix for embedding approaches
  #do a little RaggedExperiment magic
  #to allow for disjoint compartment calls
  if (is(obj, "list")) obj <- RaggedExperiment(as(obj, "GRangesList"))
  
  #condense
  fixed <- ifelse("flip.score" %in% getAssayNames(obj),
                  TRUE, FALSE)
  if (isTRUE(fixed)) se <- compactSummarizedExperiment(obj, i = "flip.score")
  if (isFALSE(fixed)) se <- compactSummarizedExperiment(obj, i = "score")
  
  #go to a matrix
  score.mat <- t(as.matrix(assay(se)))
  
  #make NAs soft zeros
  message("Making NAs soft zeros.")
  score.mat[is.na(score.mat)] <- 0
  
  #tf-idf transform
  if (isTRUE(tf.idf)) {
    message("Transforming with TF-IDF.")
    score.mat <- transformTFIDF(score.mat)
  }
  
  if (isFALSE(cluster)) return(score.mat)
  
  #dim reduction
  if (!is.null(dims) && dims <= nrow(score.mat)) {
    message("Running PCA and retaining ", dims, " dimensions.")
    pca.results <- runPCA(t(score.mat), rank = dims, BSPARAM=ExactParam())
    rownames(pca.results$rotation) <- rownames(score.mat)
  } 
  if (!is.null(dims) && dims > nrow(score.mat)) {
    stop("dims needs to be less than the number of samples/cells we have.")
    }
  
  #cluster
  message("Clustering using k=", k, " neighbors.")
  if (exists("pca.results")) {
    #compute neighbors in PCA space
    knn <- RANN::nn2(pca.results$rotation, k = k)
    adj_mat <- .getAdjacencyMatrix(pca.results$rotation, knn)
  } else {
    #compute neighbors in TF-IDF or un-transformed space
    knn <- RANN::nn2(score.mat, k = k)
    adj_mat <- .getAdjacencyMatrix(score.mat, knn)
  }
  
  #convert to graph
  adjacency_graph <- igraph::graph.adjacency(adj_mat,
                                             mode = "undirected")
  adjacency_graph <- igraph::simplify(adjacency_graph)
  
  #find communities/clusters
  message("Finding communities in the graph.")
  louvain_com <- igraph::cluster_louvain(adjacency_graph)
  com <- louvain_com$membership
  names(com) <- louvain_com$names
  if (isTRUE(cluster) && !is.null(dims)) return(list(compartment.mat=score.mat,
                                                     communities=com,
                                                     reduced_dims=pca.results))
  if (isTRUE(cluster)) return(list(compartment.mat=score.mat,
                                   communities=com))
}

#helper function for RANN results to adjacency matrix
#modeled from #https://divingintogeneticsandgenomics.rbind.io/post/clustering-scatacseq-data-the-tf-idf-way/
.getAdjacencyMatrix <- function(obj, knn) {
  #knn input is KNN results
  knn_map <- knn$nn.idx
  #create an empty matrix
  adjacency_mat <- matrix(0, nrow(obj), nrow(obj))
  #fill it in
  rownames(adjacency_mat) <- colnames(adjacency_mat) <- rownames(obj)
  for(i in seq_len(nrow(obj))) {
    adjacency_mat[i,rownames(obj)[knn_map[i,]]] <- 1
  }
  return(adjacency_mat)
}
