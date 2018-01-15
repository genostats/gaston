\name{qqplot.pvalues}
\alias{qqplot.pvalues}

\title{ QQ plot of p-values }

\description{ Draw a QQ plot of p-values }
\usage{ qqplot.pvalues(p, pch = ".", cex = 2, 
                       xlab = expression(paste("expected ", -log[10](p))),
                       ylab = expression(paste("observed ", -log[10](p))), 
                       main = "QQ plot of p-values", col.abline = "red", 
                       CB = TRUE, col.CB = gray(.9), CB.level = 0.95, ...) }

\arguments{
  \item{p}{ A vector of p-values }
  \item{pch}{ The plotting symbol }
  \item{cex}{ Size of points }
  \item{xlab}{ X-axis label }
  \item{ylab}{ Y-axis label }
  \item{main}{ Main title for the plot }
  \item{col.abline}{ Color of the line of slope 1. Set to \code{NA} to suppress. }
  \item{CB}{ \code{Logical}. If \code{TRUE}, a confidence band is included in the plot. }
  \item{col.CB}{ The color of the confidence band. }
  \item{CB.level}{ The level of the confidence band. }
  \item{...}{ Additional parameters to be passed to \code{plot} and \code{points} }
}

\details{ The QQ plot is on the \eqn{-\log_{10}}{-log10} scale, as is usual when reporting
GWAS results. The confidence band is not a global confidence region: it is the juxtaposition
of confidence intervals for each quantile.}

\seealso{ \code{\link{qqplot}}, \code{\link{plot.default}}, \code{\link{points.default}} }

\examples{
# a vector of uniform p-values
p <- runif(1e5)
qqplot.pvalues(p)
}

