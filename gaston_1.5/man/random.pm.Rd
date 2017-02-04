\name{random.pm}
\alias{random.pm}

\title{ Random square definite positive matrix }

\description{ Generate a random definite positive matrix with specified dimension }
\usage{ random.pm(n, values = n*diff(seq(0,1,length=n+1)**0.75)) }

\arguments{
  \item{n}{ Dimension of matrix }
  \item{values}{ Eigenvalues of matrix }
}

\details{
  The default for \code{values} is chosen so that the eigenvalues of the 
  resulting matrix are roughly
  similar to the eigenvalues typically observed on Genetic Relationship Matrices.
}

\value{
  A named list with members:
  \item{K}{ A \code{n x n} symmetric positive matrix }
  \item{eigen}{ The eigen decomposition of \code{K} as \code{eigen(K)} would output it }
}

\author{ Hervé Perdry and Claire Dandine-Roulland }

\seealso{ \code{\link{lmm.simu}}, \code{\link[base:eigen]{eigen}}}

\examples{
# generate a random positive matrix 
set.seed(1)
R <- random.pm(503)
str(R)
}
