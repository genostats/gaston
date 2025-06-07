#' Random square definite positive matrix
#' 
#' @description  Generate a random definite positive matrix with specified dimension
#' 
#' @param n  Dimension of matrix
#' @param values  (Optional) A numeric vector of dimension n : the eigenvalues of the matrix
#' 
#' @details
#' If \code{values} isn't given, it is chosen (deterministically)
#' so that the eigenvalues of the resulting matrix are
#' similar to eigenvalues observed on Genetic Relationship Matrices.
#' 
#' The random matrix is generated as \eqn{ U diag( values ) U' }{U \%*\% diag(values) \%*\% t(U)}
#' with \eqn{U} a random orthogonal matrix.
#' 
#' @return
#' A named list with members:
#' \item{K}{ A \code{n x n} symmetric positive matrix }
#' \item{eigen}{ The eigen decomposition of \code{K} as \code{eigen(K)} would output it }
#' @seealso  \code{\link{lmm.simu}}, \code{\link[base:eigen]{eigen}}
#' 
#' @examples
#' 
#' # generate a random positive matrix 
#' set.seed(1)
#' R <- random.pm(500)
#' str(R)
#' 
#' @export random.pm
random.pm <- function(n, values) {
  if(missing(values)) values <- grm.values(n)
  if(length(values) != n) stop("values should be of length n")
  Q <- random.ortho(n)
  K <- Q %*% (values * t(Q))
  return(list(K = K, eigen = list(values = values, vectors = Q)))
}

grm.values <- function(n) {
  kappa <- 1.63e-5
  gamma <- 0.39045 + n/13580
  delta <- 1/(1.1074 + n/175000)
  a <- 1/(1+exp( (gamma - qnorm(ppoints(n)))/delta ))
  mu.a <- mean(a)
  sd.a <- sd(a)
  rev( 1 + sqrt(kappa*n)*(a - mu.a)/sd.a )
}

random.ortho <- function(n)
  return(.Call(`_gaston_random_ortho`, PACKAGE = "gaston", n))

