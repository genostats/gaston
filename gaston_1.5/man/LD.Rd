\name{LD}
\alias{LD}

\title{ Linkage Disequilibrium }

\description{ Compute Linkage Disequilibrium (LD) between given SNPs. }
\usage{ LD(x, lim, lim2, measure = c("r2", "r", "D"), trim = TRUE) }

\arguments{
  \item{x}{ A \code{\link{bed.matrix}} }
  \item{lim}{ Range of SNPs for which the LD is computed }
  \item{lim2}{ (Optional) Second range of SNPs (see Details) }
  \item{measure}{ The LD measure }
  \item{trim}{ \code{Logical}. If \code{TRUE}, the values above 1 or below -1 are replaced by 1 and -1 respectively. }
}

\details{
  If \code{lim2} is missing, the LD is computed between all SNPs with indices between \code{lim[1]} and \code{lim[2]}; 
  else, the LD is computed between the SNPs in the range given by \code{lim} and those in the range given by \code{lim2}.
 
  Note that the LD estimates are moment estimates (which are less precise than Maximum Likelihood Estimates). 
  If \code{standardize(x) = "none"}, \code{x} will be standardized
  using \code{x@mu} and \code{x@sigma}. If \code{standardize(x) = "p"}, the moment estimates can produce \eqn{r}
  values outside of the range \eqn{[-1;1]}, hence the parameter \code{trim}. We recommend to set
  \code{standardize(x) <- "mu"} (trimming can still be necessary due to rounding errors).
}

\value{
  A matrix of LD values.
}

\author{ Hervé Perdry and Claire Dandine-Roulland }

\seealso{ \code{\link{LD.thin}},  \code{\link{LD.plot}}  }

\examples{
# Load data
data(AGT)
x <- as.bed.matrix(AGT.gen, AGT.fam, AGT.bim)

# Compute LD
ld.x <- LD(x, c(1,ncol(x)))

# Plot a tiny part of the LD matrix
LD.plot( ld.x[1:20,1:20], snp.positions = x@snps$pos[1:20] )
}

\keyword{ Linkage Disequilibrium }
