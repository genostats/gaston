\name{manhattan}
\alias{manhattan}

\title{ Manhattan plot }

\description{ Draws a Manhattan plot }
\usage{ manhattan(x, bty = "n", chrom.col = c("black", "gray50"), thinning = TRUE, ... ) }

\arguments{
  \item{x}{ A data.frame with columns named \code{chr}, \code{pos} and \code{p}. }
  \item{bty}{ Type of box to draw about the plot. Default is to draw none. }
  \item{thinning}{ \code{Logical}. If \code{TRUE}, not all points are displayed. }
  \item{chrom.col}{ Alternating colors for chromosomes. }
  \item{...}{ Graphical parameters to be passed to \code{plot}. }
}

\details{If there is only one chromosome value in \code{x$chr}, the x-axis will be labeled with the SNP
position. In the general case, the x-axis is labeled with the chromosome name and the color
of the points alternates between the colors in \code{chrom.col}.

The default value \code{bty = "n"} should give the best result for GWAS Manhattan plots. 
See \code{\link{par}} for other possible values of \code{bty} and their meaning.

The thinning procedure suppress some points to avoid generating too heavy graphs. The user 
should check that setting \code{thinning = FALSE} does not change the final aspect of the 
plot. 
}

\value{An invisible copy of \code{x} is returned, in which a column \code{coord} has been added  
if there is more than one chromosome value in \code{x$chr}. This column contains the x-coordinates of each
SNP on the plot, and should prove helpful to annotate it. }

\seealso{ \code{\link{association.test}}, \code{\link{qqplot.pvalues}}, 
          \code{\link{par}}, \code{\link{plot.default}}, \code{\link{points.default}} }

