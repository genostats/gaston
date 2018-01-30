\name{qqplot.pvalues}
\alias{qqplot.pvalues}

\title{ QQ plot of p-values }

\description{ Draws a QQ plot of p-values }
\usage{ qqplot.pvalues(p, col.abline = "red", CB = TRUE, col.CB = "gray80", 
                       CB.level = 0.95, thinning = TRUE, ...) }

\arguments{
  \item{p}{ A vector of p-values, or a data.frame with a column named \code{p} }
  \item{col.abline}{ Color of the line of slope 1. Set to \code{NA} to suppress. }
  \item{CB}{ \code{Logical}. If \code{TRUE}, a confidence band is included in the plot. }
  \item{col.CB}{ The color of the confidence band. }
  \item{CB.level}{ The level of the confidence band. }
  \item{thinning}{ \code{Logical}. If \code{TRUE}, not all points are displayed. }
  \item{...}{ Graphical parameters to be passed to \code{plot} and \code{points} }
}

\details{ The QQ plot is on the \eqn{-\log_{10}}{-log10} scale, as is usual when reporting
GWAS results. 

The confidence band is not a global confidence region: it is the mere juxtaposition 
of confidence intervals for each quantile. Moreover it assumes independance of the
p-values, an hypothesis hich is false for the p-values resulting from an association 
test in presence of linkage disequilibrium. Therefore, the probability that some of the
points lie outsite of this band is greater that \code{CB.level}.

The thinning procedure suppress some points to avoid generating too heavy graphs. The user 
should check that setting \code{thinning = FALSE} does not change the final aspect of the 
QQ plot. 
}

\seealso{ \code{\link{association.test}}, \code{\link{manhattan}}, \code{\link{qqplot}}, 
          \code{\link{plot.default}}, \code{\link{points.default}} }

\examples{
# a vector of uniform p-values
p <- runif(1e6)
qqplot.pvalues(p)
# if we don't thin the points, using pch = "." is advised
qqplot.pvalues(p, pch = ".", cex = 2, thinning = FALSE)
}

