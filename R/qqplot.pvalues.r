#' QQ plot of p-values
#' 
#' @description  Draws a QQ plot of p-values
#' 
#' @param p  A vector of p-values, or a data.frame with a column named \code{p}
#' @param col.abline  Color of the line of slope 1. Set to \code{NA} to suppress.
#' @param CB  \code{Logical}. If \code{TRUE}, a confidence band is included in the plot.
#' @param col.CB  The color of the confidence band.
#' @param CB.level  The level of the confidence band.
#' @param thinning  \code{Logical}. If \code{TRUE}, not all points are displayed.
#' @param ...  Graphical parameters to be passed to \code{plot} and \code{points}
#' 
#' 
#' @details  The QQ plot is on the \eqn{-\log_{10}}{-log10} scale, as is usual when reporting
#' GWAS results.
#' 
#' The confidence band is not a global confidence region: it is the mere juxtaposition
#' of confidence intervals for each quantile. Moreover it assumes independance of the
#' p-values, an hypothesis hich is false for the p-values resulting from an association
#' test in presence of linkage disequilibrium. Therefore, the probability that some of the
#' points lie outsite of this band is greater that \code{CB.level}.
#' 
#' The thinning procedure suppress some points to avoid generating too heavy graphs. The user
#' should check that setting \code{thinning = FALSE} does not change the final aspect of the
#' QQ plot.
#' 
#' 
#' 
#' @seealso  \code{\link{association.test}}, \code{\link{manhattan}}, \code{\link{qqplot}},
#' \code{\link{plot.default}}, \code{\link{points.default}}
#' 
#' 
#' @examples
#' 
#' # a vector of uniform p-values
#' p <- runif(1e6)
#' qqplot.pvalues(p)
#' # if we don't thin the points, using pch = "." is advised
#' qqplot.pvalues(p, pch = ".", cex = 2, thinning = FALSE)
#' 
#' @export qqplot.pvalues
qqplot.pvalues <- function(p, col.abline = "red", CB = TRUE, col.CB = "gray80", CB.level = 0.95,
                           thinning = TRUE, ...) {
  if(is.list(p)) { # ok pour les data frame aussi
    if(is.null(p$p))
      stop("No p-values were found")
    p <- p$p
  }
  p <- p[ !is.na(p) ] # on supprime ces valeurs silencieusement...
  
  if(any(p>1) | any(p<0))
    stop("p-values should be in [0,1]\n")

  w <- (p == 0)
  if(any(w)) { # mais ça on avertit
    warning("There are ", sum(w), " zero p-values that won't be displayed")
    p <- p[!w]
  }

  args <- list(...)
  if(is.null(args$xlab))
    args$xlab <- expression(paste("expected ", -log[10](p)))
  if(is.null(args$ylab))
    args$ylab <- expression(paste("observed ", -log[10](p)))
  if(is.null(args$main))
    args$main <- "QQ plot of p-values"

  n <- length(p)
  args$type <- "n"
  args$x <- -log10( c(1,n) / (n+1) )
  args$y <- -log10( range(p) )
  do.call( plot, args )

  # confidence interval
  k <- n * 10**seq(-log10(n), 0, length = 200)
  hi <- -log10(qbeta( 0.5-CB.level/2, k, n+1-k))
  lo <- -log10(qbeta( 0.5+CB.level/2, k, n+1-k))
  polygon( -log10(c(k, rev(k))/(n+1)), c(lo, rev(hi)), col = col.CB, border= col.CB )

  # diag line with real limits
  segments(0, 0, -log10(1/n), -log10(1/n), col = col.abline) 

  # the points
  expected <- -log10( (n:1)/(n+1) )
  observed <- sort(-log10( p ))

  if(thinning) {
    w <- manhattan.thinning(expected, observed, 10000, 10000)
    args$x <- expected[w]
    args$y <- observed[w] 
  } else {
    args$x <- expected
    args$y <- observed
  }
  args$type <- "p"
  do.call( points, args )
}

