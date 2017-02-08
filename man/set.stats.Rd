\name{set.stats}
\alias{set.stats}
\alias{set.stats.ped}
\alias{set.stats.snps}

\title{ Basic statistics for a \code{\link{bed.matrix}} }

\description{
  Return an updated \code{\link{bed.matrix}} with new variables for 
  several basic statistics.
}
\usage{
set.stats(x, set.p = TRUE, set.mu_sigma = TRUE, 
          verbose = getOption("gaston.verbose",TRUE))

set.stats.snps(x, set.p = TRUE, set.mu_sigma = TRUE, 
               verbose = getOption("gaston.verbose",TRUE)) 

set.stats.ped(x, verbose = getOption("gaston.verbose",TRUE))
}

\arguments{
  \item{x}{ A \code{\link{bed.matrix}} }
  \item{set.p}{ If \code{TRUE}, \code{x@p} is updated }
  \item{set.mu_sigma}{ If \code{TRUE}, \code{x@mu} and \code{x@sigma} are updated }
 \item{verbose}{ If \code{TRUE}, display information on the function actions }
}

\details{
  \code{set.stats} is called by default by all functions that create a bed.matrix, unless
  the global option \code{gaston.auto.set.stats} is \code{FALSE} (cf example below).

  \code{set.stats} and \code{set.stats.ped} update \code{x@ped}, adding the following variables: 
\itemize{
  \item \code{N0}, \code{N1}, \code{N2} and \code{NAs} give for each individual the number of 
        autosomal SNPs with a genotype equal to 0, 1, 2 and missing, respectively
  \item \code{N0.x}, \code{N1.x}, \code{N2.x} and \code{NAs.x} idem for chromosome X
  \item \code{N0.y}, \code{N1.y}, \code{N2.y} and \code{NAs.y} idem for chromosome Y
  \item \code{N0.mt}, \code{N1.mt}, \code{N2.mt} and \code{NAs.mt} idem for mitochondrial SNPs
  \item \code{callrate}, \code{callrate.x}, \code{callrate.y}, \code{callrate.mt} 
        is the individual callrate for autosomal, X, Y, mitochondrial SNPs
  \item \code{hz}, \code{hz.x}, \code{hz.y}, \code{hz.mt} is the individual heterozygosity 
        for autosomal, X, Y, mitochondrial SNPs
  }

  \code{set.stats} and \code{set.stats.snps} update \code{x@snps}, adding the following variables:
\itemize{
  \item \code{N0}, \code{N1}, \code{N2} and \code{NAs} give for each SNP the number of individuals
  with a genotype equal to 0, 1, 2 and missing, respectively
  \item \code{N0.f}, \code{N1.f}, \code{N2.f} and \code{NAs.f} give, only for SNPs on chromosome X, 
  the number of female individuals with a genotype equal to 0, 1, 2 and missing, respectively
  \item \code{callrate} is the SNP callrate (for Y linked SNPs, the callrate is computed usin
     males only).
  \item \code{maf} is the Minor Allele Frequency 
  \item \code{hz} is the SNP heterozygosity (for X linked SNPs, the heterozygosity is computed
     using females only).
  }

  If \code{set.p = TRUE}, \code{x@p} is updated with the alternate allele frequency.

  If \code{set.mu_sigma = TRUE}, \code{x@mu} is updated with the genotype mean (equal to \code{2*x@p})
  and \code{x@sigma} with the genotype standard deviation (should be approximately \code{sqrt(2*x@p*(1-x@p))}
  under Hardy-Weinberg Equilibrium).
}

\value{ 
  A \code{\link{bed.matrix}} similar to \code{x}, with slots updated as described above.
}

\author{ Hervé Perdry and Claire Dandine-Roulland }

\seealso{ \code{\link{set.hwe}}, \code{\link{set.genomic.sex}} }

\examples{
# Disable auto set stats :
options(gaston.auto.set.stats = FALSE)

# Load data
data(TTN)
x <- as.bed.matrix(TTN.gen, TTN.fam, TTN.bim)
str(x@ped)
str(x@snps)

# Compute statistics
x <- set.stats(x)
str(x@ped)
str(x@snps)

# restore default behavior
options(gaston.auto.set.stats = TRUE)
}

\keyword{ Genetic }
\keyword{ Statistics }
