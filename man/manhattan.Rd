\name{manhattan}
\alias{manhattan}

\title{ Manhattan plot }

\description{ Draws a Manhattan plot }
\usage{ manhattan(x, bty = "n", thinning = TRUE, thinning.step = 1e-5, ... ) }

\arguments{
  \item{x}{ A data.frame with columns named \code{chr}, \code{pos} and \code{p}. }
  \item{bty}{ Type of box to draw about the plot. See \code{\link{par}} for possible values. } 
  \item{thinning}{ \code{Logical}. If \code{TRUE}, not all points are displayed. }
  \item{thinning.step}{ The minimal distance between thinned points (measured on the y-axis). }
  \item{...}{ Graphical parameters to be passed to \code{plot}. }
}

\details{If there is only one chromosome value in \code{x$chr}, the x-axis will be labeled with the SNP
position. In the general case, the x-axis will be labeled with the chromosome name.

The thinning procedure suppress some points to avoid generating too heavy graphs. The user 
should check that setting \code{thinning = FALSE} does not change the final aspect of the 
QQ plot. Smaller values of \code{thinning.step} result in more points plotted. For publication
quality graphs we recommend \code{thinning = FALSE}.
}

\value{An invisible copy of \code{x} is returned, in which a column \code{coord} has been added  
if there is more than one chromosome value in \code{x$chr}. This column contains the x-coordinates of each
SNP, and should prove helpful to annotate the graph. }

\seealso{ \code{\link{association.test}}, \code{\link{qqplot.pvalues}}, 
          \code{\link{par}}, \code{\link{plot.default}}, \code{\link{points.default}} }

