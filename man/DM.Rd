\name{DM}
\alias{DM}

\title{Dominance Matrix}

\description{
  Compute the Dominance Matrix}

\usage{ DM(x, which.snps, autosome.only = TRUE, chunk = 1L) }

\arguments{
  \item{x}{ A \code{\link{bed.matrix}} }
  \item{which.snps}{ Logical vector, giving which snps to use in the computation.  The default is to use all autosomal SNPs}
  \item{autosome.only}{ If \code{TRUE}, only autosomal SNPs will be considered. }
  \item{chunk}{ Parameter for the parallelization: how many SNPs are treated by each task}
}

\details{
  The Dominance Matrix (DM) gives for each pair of individuals an estimation of 
  their probability of sharing two alleles Identical By Descent. 

  It is computed by a moment estimator,
  \eqn{{1\over q} ZZ'}{ZZ'/q} with \eqn{Z} the matrix with entries
  \eqn{p \over 1-p}{p/(1-p)}, \eqn{-1}, \eqn{1-p \over p}{(1-p)/p} according to the
  values 0, 1, 2 in the genotype matrix \code{x} (here \eqn{p} is the
  frequency of the alternate allele, and \eqn{q} is the number of SNPs
  (\code{ncol(x)}).
}

\value{ A symmetric square matrix of dimension equal to the number of individuals.
Each entry can be interpreted as an estimated probability of sharing two alleles IBD 
— as it is a moment estimator, the value can (and will) fall outside of the range (0,1). }

\seealso{ \code{\link{GRM}}, \code{\link{reshape.GRM}} }
\examples{
# load chr2 data set (~10k SNPs in low LD)
x <- read.bed.matrix( system.file("extdata", "chr2.bed", package="gaston") )

# Compute Dominance Matrix
D <- DM(x)
dim(D)

}

\keyword{ Dominance Matrix }
