\name{LD.thin}
\alias{LD.thin}

\title{ LD thinning }

\description{ Select SNPs in LD below a given threshold. }
\usage{
LD.thin(x, threshold, max.dist = 250e3, beg = 1, end = ncol(x),
        which.snps, dist.unit = c("bases", "indices", "cM"), 
        extract = TRUE, keep = c("left", "right", "random"))
}

\arguments{
  \item{x}{ A \code{\link{bed.matrix}} }
  \item{threshold}{ The maximum LD (measured by \eqn{r^2}) between SNPs }
  \item{max.dist}{ The maximum distance for which the LD is computed }
  \item{beg}{ The index of the first SNP to consider }
  \item{end}{ The index of the last SNP to consider }
  \item{which.snps}{ Logical vector, giving which SNPs are considerd. The default is to use all SNPs}
  \item{dist.unit}{ Distance unit in \code{max.dist} }
  \item{extract}{ A \code{logical} indicating whether the function return a \code{bed.matrix} (\code{TRUE}) 
                  or a logical vector indicating which SNPs are selected (\code{FALSE}) }
  \item{keep}{ Which SNP is selected in a pair with LD above \code{threshold} }
}

\details{
  The SNPs to keep are selected by a greedy algorithm. The LD is computed only for SNP pairs for which distance is inferior to 
  \code{max.dist}, expressed in number of bases if \code{dist.unit = "bases"}, in number of SNPs if \code{dist.unit = "indices"},
  or in centiMorgan if \code{dist.unit = "cM"}. 
 
  The argument \code{which.snps} allows to consider only a subset of SNPs.

  The algorithm tries to keep the largest possible number of SNPs: it is not appropriate to select tag-SNPs.
}

\value{
  If \code{extract = TRUE}, a \code{\link{bed.matrix}} extracted from \code{x} with SNPs in pairwise LD below the given threshold.
  If \code{extract = FALSE}, a logical vector of length \code{end - beg + 1}, where \code{TRUE} indicates that 
  the corresponding SNPs is selected.
}

\author{ Hervé Perdry and Claire Dandine-Roulland }

\seealso{ \code{\link{LD}}, \code{\link{set.dist}} }

\examples{
# Load data
data(TTN)
x <- as.bed.matrix(TTN.gen, TTN.fam, TTN.bim)

# Select SNPs in LD r^2 < 0.4, max.dist = 500 kb
y <- LD.thin(x, threshold = 0.4, max.dist = 500e3)
y

# Verifies that there is no SNP pair with LD r^2 > 0.4
# (note that the matrix ld.y has ones on the diagonal)
ld.y <- LD( y, lim = c(1, ncol(y)) )
sum( ld.y > 0.4 )  
}

\keyword{ Linkage Disequilibrium }
