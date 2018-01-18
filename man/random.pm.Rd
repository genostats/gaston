\name{random.pm}
\alias{random.pm}

\title{ Random square definite positive matrix }

\description{ Generate a random definite positive matrix with specified dimension }
\usage{ random.pm(n, values) }

\arguments{
  \item{n}{ Dimension of matrix }
  \item{values}{ (Optional) A numeric vector of dimension n : the eigenvalues of the matrix }
}

\details{
  If \code{values} isn't given, it is chosen (deterministically) 
  so that the eigenvalues of the resulting matrix are 
  similar to eigenvalues observed on Genetic Relationship Matrices.

  The random matrix is generated as \eqn{ U diag( values ) U' }{U \%*\% diag(values) \%*\% t(U)}
  with \eqn{U} a random orthogonal matrix.
}

\value{
  A named list with members:
  \item{K}{ A \code{n x n} symmetric positive matrix }
  \item{eigen}{ The eigen decomposition of \code{K} as \code{eigen(K)} would output it }
}

\seealso{ \code{\link{lmm.simu}}, \code{\link[base:eigen]{eigen}}}

\examples{
# generate a random positive matrix 
set.seed(1)
R <- random.pm(500)
str(R)
}
