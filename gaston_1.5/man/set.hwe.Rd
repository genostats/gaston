\name{set.hwe}
\alias{set.hwe}

\title{ Hardy-Weinberg Equilibrium }

\description{
  Returns an updated \code{\link{bed.matrix}} with a new variable for the \eqn{p}-values of an
  Hardy-Weinberg Equilibrium test.
}
\usage{ set.hwe(x, method = c("chisquare", "exact"),
           verbose = getOption("gaston.verbose", TRUE))
}

\arguments{
  \item{x}{ A \code{\link{bed.matrix}} }
  \item{method}{ The method to use, either "chisquare" or "exact"  }
 \item{verbose}{ If \code{TRUE}, display information on the function actions }
}

\details{
Two tests of Hardy-Weinberg Equilibrium are proposed:
\itemize{
  \item if \code{method = "chisquare"}, the good old Chi-square test
  \item if \code{method = "exact"}, Haldane's exact test (see Wigginton et al)
  }

The function \code{set.stats} will be called first if necessary.

The \eqn{p}-value is set to \eqn{1.0} for SNPs on chromosomes Y and MT. For SNPs on
chromosomes X, currently, the test is performed using only the genotypic counts of women.
}

\value{
  A \code{\link{bed.matrix}} similar to \code{x}, with a new variable \code{x@snps$hwe}
  containing the \eqn{p}-values for each SNP.
}

\references{ Wigginton, J. E., Cutler, D. J., & Abecasis, G. R. (2005). \emph{A note on exact tests of Hardy-Weinberg equilibrium}. The American Journal of Human Genetics, \bold{76(5)}, \bold{887-893} }

\author{Hervé Perdry and Claire Dandine-Roulland}

\seealso{ \code{\link{set.stats}}, \code{\link{set.genomic.sex}} }

\examples{
# Load data
data(LCT)
x <- as.bed.matrix(LCT.gen, LCT.fam, LCT.bim)

# Compute Hardy-Weinberg p-values
x <- set.hwe(x)
head( x@snps[,c("id","hwe")] )
}
\keyword{ Hardy-Weinberg }
\keyword{ P-value }
