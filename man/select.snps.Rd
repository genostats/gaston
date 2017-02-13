\name{select.snps}
\alias{select.snps}

\title{ Subsetting from a \code{\link{bed.matrix}} }

\description{
Returns subset of SNPs satisfying a condition.
}

\usage{ select.snps(x, condition) }

\arguments{
  \item{x}{ A \code{\link{bed.matrix}} }
  \item{condition}{ Condition used to select SNPs }
}

\details{
The conditions can involve global variables and all variables defined
in the data frame \code{x@snps}, in particular
\itemize{
  \item \code{chr}, \code{id}, \code{dist}, \code{pos}, \code{A1}, \code{A2}
  \item If basic stats have been computed (see \code{\link{set.stats}}), \code{N0}, \code{N1}, \code{N2}, \code{NAs}, \code{callrate}, \code{maf}, \code{hz}, etc.
  \item If Hardy-Weinberg Equilibrium test has been performed (see \code{\link{set.hwe}}), \code{hwe}.
  }

If some condition evaluate to \code{NA} (e.g. \code{maf > 0} when \code{maf} is undefined for some SNPs), 
a warning is issued and the corresponding SNPs are removed.
}

\value{ A \code{\link{bed.matrix}} similar to \code{x}, containing the selected SNPs only }

\author{ Hervé Perdry and Claire Dandine-Roulland }

\seealso{ \code{\link{select.snps}}, \code{\link{set.stats}}, \code{\link{set.hwe}} }

\examples{
# Load data
data(LCT)
x <- as.bed.matrix(LCT.gen, LCT.fam, LCT.bim)

# Select SNPs with a maf > 5%
y <- select.snps(x, maf > 0.05)
y
}
