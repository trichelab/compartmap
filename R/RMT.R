#' Random matrix theory plot, stolen from `covmat`, which fell off of CRAN
#' 
#' Eigenvalue plot
#' 
#' @details
#' Plots eigenvalues of the correlation matrix and overlays the Marchenko-
#' Pastur density on top of it. There is a sharp cutoff for the density. We
#' are concerned with eigenvalues beyond this cutoff. Parameters used for 
#' plotting are added to the plot. 
#' 
#' TODO: replace the example with one from compartmap
#' TODO: consider making RMT an S4 class someday 
#' 
#' @param x     model of type RMT obtained by fitting an RMT model to the data
#' @param y     unused
#' @param ...   unused
#'
#' @author Rohit Arora
#' 
#' @examples 
#' \dontrun{
#'  data("largereturn")
#'  model <- estRMT(largesymdata)
#'  plot(model)
#' }
#' 
#' @method plot RMT
#' 
#' @import timeSeries
#' @import ggplot2
#' @import scales
#' 
#' @export
plot.RMT <- function(x, y, ...){
    
  lambdas <- x$eigVals
  Q <- x$Q
  sigma.sq <- x$var 
  lambda.max <- x$lambdascutoff 
    
  p <- ggplot(data=data.frame(lambdas)) + 
         geom_histogram(aes_string(x=lambdas, 
                                   y='..density..'),
                        breaks=seq(min(lambdas)-1, 1+max(lambdas),0.5), 
                        colour="black", 
                        fill="white") +
         stat_function(fun=dmp, 
                       args=list(svr=Q, var=sigma.sq), 
                       aes(colour='MP density')) + 
         xlab("Eigenvalues") +
         labs(title="Actual vs Fitted Marchenko-Pastur") + 
         ylim(0,1.5) + 
         theme(plot.title=element_text(size=20, 
                                       face="bold", 
                                       vjust=1),
               axis.title=element_text(size=14,
                                       face="bold")) + 
         annotate('text', 
                  x=10, 
                  y=0.9, 
                  label=paste("sigma^{2} == ", round(sigma.sq,3)), 
                  parse=TRUE) +
         annotate('text', 
                  x=10, 
                  y=1, 
                  label=paste("Q == ", round(Q,3)), 
                  parse=TRUE) + 
         annotate('text', x=10, y=0.78, 
                  label=paste("lambda[max] ==", round(lambda.max,3)), 
                  parse=TRUE) + 
         scale_colour_manual("", values=c("red")) + 
         NULL # standard practice for long ggplots
    
    options(warn=-1)
    print(p) # wtf
    options(warn=0)
    p
}
