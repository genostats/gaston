\name{GRM}
\alias{GRM}

\title{ Genetic Relationship Matrix}

\description{
  Compute the Genetic Relationship Matrix}

\usage{ GRM(x, which.snps, autosome.only = TRUE, chunk = 1L) }

\arguments{
  \item{x}{ A \code{\link{bed.matrix}} }
  \item{which.snps}{ Logical vector, giving which snps to use in the computation.  The default is to use all autosomal SNPs}
  \item{autosome.only}{ If \code{TRUE}, only autosomal SNPs will be considered. }
  \item{chunk}{ Parameter for the parallelization: how many SNPs are treated by each task}
}

\details{
  The Genetic Relationship Matrix (GRM) is computed by the formula \eqn{{1\over q}XX'}{XX'/q}, 
  with \eqn{X} the standardized genotype matrix and \eqn{q} the number of SNPs
  (\code{ncol(x)}).

  If \code{x} is not standardized before this computation, the function
  will use \code{standardize(x) <- "p"} by default.
}

\value{ The GRM is a symmetric square matrix of dimension equal to the number of individuals.
Each entry can be interpreted as an estimated kinship coefficient between individuals, although some 
authors might disagree. Note in particular that some entries will be negative.}

\author{ Hervé Perdry and Claire Dandine-Roulland }

\seealso{ \code{\link{DM}}, \code{\link{reshape.GRM}}, \code{\link{lmm.aireml}}, \code{\link{lmm.diago}}, \code{\link{standardize}}, \code{\link{bed.loadings}} }
\examples{
# load chr2 data set (~10k SNPs in low LD)
x <- read.bed.matrix( system.file("extdata", "chr2.bed", package="gaston") )

# Compute Genetic Relationship Matrix
K <- GRM(x)
dim(K)
}

\keyword{ Genetic Relationship Matrix }
