#' Contour plot for two parameters likelihood
#' 
#' @description  Create a contour plot (superimposed with a heat map)
#' 
#' @param x,y,z  As in \code{contour}
#' @param levels  As in \code{contour}. If \code{NULL}, the function computes appropriate levels.
#' @param nlevels  As in \code{contour}
#' @param heat  If \code{TRUE}, a heat map is superimposed to the contour plot
#' @param col.heat  Vector of heat colors
#' @param \dots  Additional arguments to \code{image} and \code{contour}
#' 
#' @details  This function is a wrapper for \code{contour}, with a different method to compute
#' a default value for levels. If \code{heat = TRUE}, a heatmap produced by \code{image} is added to the plot.
#' See \code{\link{contour}} for details on parameters.
#' 
#' @seealso  \code{\link{lmm.diago.likelihood}}, \code{\link[graphics:contour]{contour}}, \code{\link[graphics:image]{image}}
#' 
#' @keywords  Heat map
#' @examples
#' 
#' data(AGT)
#' x <- as.bed.matrix(AGT.gen, AGT.fam, AGT.bim)
#' 
#' # Compute Genetic Relationship Matrix
#' K <- GRM(x)
#' 
#' # eigen decomposition of K
#' eiK <- eigen(K)
#' 
#' # simulate a phenotype
#' set.seed(1)
#' y <- 1 + lmm.simu(tau = 1, sigma2 = 2, eigenK = eiK)$y
#' 
#' # Likelihood
#' TAU <- seq(0.5,2.5,length=30)
#' S2 <- seq(1,3,length=30)
#' lik1 <- lmm.diago.likelihood(tau = TAU, s2 = S2, Y = y, eigenK = eiK)
#' lik.contour(TAU, S2, lik1, heat = TRUE, xlab = "tau", ylab = "sigma^2")
#' 
#' @export lik.contour
lik.contour <- function(x, y, z, levels=NULL, nlevels=11, heat = TRUE, col.heat=NULL, ...)
{
  if(is.null(levels))
  {
    zz <- z[z < Inf & z > -Inf & !is.na(z)]
    levels <- quantile( zz , c(0.001,0.007, seq( 0.05, 0.95, length=nlevels-4),0.993, 0.999) );
    d <- max(-log(min(abs(diff(levels))))/log(10),0);
    levels <- round(levels,d+1);
  }
  if(is.null(col.heat)) col.heat <- heat.colors(length(levels)-1)
  if(heat) image(x,y,z, breaks= sort(levels), col=col.heat,...)
  contour(x,y,z, levels=levels,add=heat,...)
  invisible( list(x=x, y=y, z=z) )
}
