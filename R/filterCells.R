#' Filter cells using either k-means or Dirichlet process means clustering of sparsity metrics
#'
#' @name filterCells
#'
#' @param sparsity.mat A matrix of summarized sparsity measures
#' @param rse.obj The unfiltered RangedSummarizedExperiment object
#' @param cluster.method Clustering method to use (default: kmeans)
#' @param clusters How many clusters to generate
#' @param plot.data Whether to plot the data
#'
#' @return A filtered RangedSummarizedExperiment object and/or plot of the filtered data
#' @export
#' @import RColorBrewer
#' @import grid
#' @import SummarizedExperiment
#'
#' @examples
#'

filterCells <- function(sparsity.mat, rse.obj, cluster.method = c("kmeans", "dpmeans"),
                        clusters = 2, plot.data = FALSE) {
  #this function is to filter out poor/low quality cells given the sparsity metrics
  #the underlying idea is that we will have 2 clusters of cells
  #1. cells with sufficient signal to use in our inference
  #2. cells with low/too sparse of signal genome-wide
  #input comes from getBinSparsity()
  #output will be a filtered RangeSummarizedExperiment object
  #ideally we can do this with lambda means - but not implemented yet...
  
  #check if the input is a per locus matrix and summarize if needed
  # if (is(sparsity.mat, "list")) {
  #   if (!("scaled.sparsity.mat" %in% names(sparsity.mat[[1]]))) {
  #     message("Summarizing sparsity measures.")
  #     if ("sparsity.mat" %in% names(sparsity.mat)) {
  #       sparsity.mat$sparsity.mat <- colSums(sparsity.mat$sparsity.mat) / nrow(sparsity.mat$sparsity.mat)
  #       single.chr <- TRUE
  #     } else {
  #       sparsity.mat <- lapply(sparsity.mat, function(x) {
  #         if ("sparsity.mat" %in% names(x)){
  #           scaled.sparsity.mat <- colSums(x[1]) / nrow(x[1])
  #           } else {
  #             scaled.sparsity.mat <- colSums(x) / nrow(x)
  #           }
  #       })
  #       single.chr <- FALSE
  #     }
  #   } else { single.chr <- FALSE }
  # } else { single.chr <- FALSE } 
  # 
  # #assumes that the object was done using getBinMatrix and lapply
  # if (is(sparsity.mat, "list") & single.chr) {
  #   sparsity.mat <- as.matrix(sparsity.mat[1])
  # } else if (is(sparsity.mat, "list")){
  #   sparsity.mat.tmp <- lapply(sparsity.mat, function(x) x[1])
  #   sparsity.mat <- plyr::ldply(sparsity.mat.tmp, "rbind")
  #   }
  
  #check if the rse.obj is a RangedSummarizedExperiment object
  if (!is(rse.obj, "RangedSummarizedExperiment")) stop(rse.obj, " does not look like a RangedSummarizedExperiment object.")
  
  #get the clustering method
  cluster.method <- match.arg(cluster.method)
  
  #perform the clustering
  message("Clustering the data using: ", cluster.method)
  if (cluster.method == "kmeans") cluster.results <- stats::kmeans(t(sparsity.mat), centers = clusters)
  if (cluster.method == "dpmeans") stop("DPmeans is not implemented yet.")
  
  #filter the data
  if (!is.null(colnames(rse.obj))) {
    colname.match <- match(colnames(rse.obj), names(cluster.results$cluster))
    assay(rse.obj) <- assay(rse.obj)[,colname.match]
    }
  message("Filtering cells with too sparse of signal from ", deparse(substitute(rse.obj)), ".")

  #get the cluster to filter by sorting the per-sample sparsity measures and picking the lowest one
  #and querying which cluster that cell/sample belongs to
  cluster.to.filter <- cluster.results$cluster[names(sort(rowSums(t(sparsity.mat))))[1]]
  message("Filtering out cells belonging to cluster ", as.character(cluster.to.filter))
  filter.vec <- ifelse(cluster.results$cluster == cluster.to.filter, FALSE, TRUE) #make sure logic is inverted
  filtered.data <- rse.obj[,filter.vec]
  message("Filtered out ", table(filter.vec)["FALSE"], " cells.")
  message("Kept ", table(filter.vec)["TRUE"], " cells with sufficient signal.")
  
  #plot the data
  if (plot.data) {
    message("Plotting results.")
    pca_sparsity <- prcomp(t(sparsity.mat))
    pr_comps <- data.frame(pca_sparsity$x)
    
    # Combine for plotting
    pr_comps$clusters <- factor(cluster.results$cluster)
    
    cols <- colorRampPalette(brewer.pal(3, "Paired"))
    
    pca_plot <- ggplot(pr_comps, aes(x=PC1, y=PC2, color=clusters)) + 
      geom_point(size=3.5) + 
      ylim(-40, 40) +
      xlim(-40, 40) +
      scale_color_manual(values=cols(2)) + 
      theme_bw(10)
    
    # Plot percent variation explained
    prop_var <- data.frame(t(summary(pca_sparsity)$importance))
    names(prop_var) = c('sd', 'prop', 'cum')
    prop_var$num = 1:nrow(prop_var)
    
    var_plot <- ggplot(prop_var, aes(x=num, y=prop)) + 
      geom_point(size=3.5) + 
      geom_line() + 
      scale_x_continuous(limits = c(1, 22), breaks = 1:22) +
      xlab("Principal Component") + 
      ylab("Prop. of Variance") +
      theme_bw(10) +
      theme(
        axis.title.y = element_text(vjust=1),
        plot.margin = unit(c(0,0,0,6), "mm")
      )
    
    vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
    
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(4, 100)))
    print(pca_plot, vp = vplayout(1:3, 3:100))
    print(var_plot, vp = vplayout(4, 1:92))
  }
  
  return(filtered.data)
}
