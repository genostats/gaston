\name{Tests}
\alias{test.inds}
\alias{test.snps}
\alias{which.inds}
\alias{which.snps}

\title{ Evaluation of a condition on SNPS or individuals in a \code{\link{bed.matrix}} }

\description{
Evaluate a condition and return logical vector or indices
}

\usage{ test.snps(x, condition, na.to.false = TRUE) 
        test.inds(x, condition, na.to.false = TRUE) 
        which.snps(x, condition) 
        which.inds(x, condition) }

\arguments{
  \item{x}{ A \code{\link{bed.matrix}} }
  \item{condition}{ Condition used to select SNPs }
  \item{na.to.false}{ If \code{TRUE}, \code{NA}s are replaced by \code{FALSE} }
}

\details{
The conditions can involve global variables and all variables defined
in the data frame \code{x@snps}, in particular for \code{test.snps} and \code{which.snps}
\itemize{
  \item \code{chr}, \code{id}, \code{dist}, \code{pos}, \code{A1}, \code{A2}
  \item If basic stats have been computed (see \code{\link{set.stats}}), \code{N0}, \code{N1}, \code{N2}, \code{NAs}, \code{callrate}, \code{maf}, \code{hz}, etc.
  \item If Hardy-Weinberg Equilibrium test has been performed (see \code{\link{set.hwe}}), \code{hwe}.
  }
and for  \code{test.inds} and \code{which.inds}
\itemize{
  \item \code{famid}, \code{id}, \code{father}, \code{mother}, \code{sex}, \code{pheno}
  \item If basic stats have been computed (see \code{\link{set.stats}}),
        \code{N0}, \code{N1}, \code{N2}, \code{NAs}, \code{callrate}, etc.
  }
}
\value{ \code{test.snps} and \code{test.inds} return a logical vector of length \code{ncol(x)} and \code{nrow(x)} respectively. \code{which.snps(x, condition)} is
equivalent to \code{which(test.snps(x, condition))} and \code{which.inds(x, condition)} to \code{which(test.inds(x, condition))}.
}

\seealso{ \code{\link{select.snps}}, \code{\link{select.inds}}, \code{\link{set.stats}}, \code{\link{set.hwe}} }

\examples{
# Load data
data(LCT)
x <- as.bed.matrix(LCT.gen, LCT.fam, LCT.bim)

# SNPs and individuals with a callrate < 100%
w <- test.snps(x, callrate < 1)
table(w)
which.snps(x, callrate < 1)
which.inds(x, callrate < 1)
}
