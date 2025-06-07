#' Likelihood of a linear mixed model
#' 
#' @description
#' Compute the Restricted or the Full Likelihood of a linear mixed model.
#' 
#' @param Y  Phenotype vector
#' @param X  Covariable matrix
#' @param K  A positive definite matrix or a \code{list} of such matrices
#' @param tau  Value(s) of parameter(s) \eqn{\tau}{tau}
#' @param s2   Value of parameter \eqn{\sigma^2}{s2}
#' @param h2   Value(s) of heritability
#' 
#' 
#' @details
#' 
#' Theses function respectively compute the Restricted and the Profile Likelihood under the linear
#' mixed model
#' \deqn{ Y = X\beta + \omega_1 + \ldots + \omega_k + \varepsilon }{ Y =X beta + omega_1 + ... + omega_k + epsilon }
#' with \eqn{ \omega_i \sim N(0,\tau_i K_i) }{omega_i ~ N(0, tau_i K_i)} for \eqn{ i \in 1, \dots,k } and
#' \eqn{ \varepsilon \sim N(0,\sigma^2 I_n) }{epsilon ~ N(0, sigma^2 I_n)}.
#' 
#' The variance matrices \eqn{K_1}{K_1}, ..., \eqn{K_k}{K_k}, are specified through the parameter \code{K}.
#' The parameter \code{tau} should be a vector of length \eqn{k}.
#' 
#' The function \code{lmm.restricted.likelihood} computes the restricted
#' likelihood for the given values of \eqn{\tau}{tau} and \eqn{\sigma^2}{s2}.
#' Whenever \eqn{k = 1}, it is similar to \code{lmm.diago.likelihood(tau, s2, Y = Y, X = X, eigenK = eigen(K))}
#' which should be prefered (with a preliminary computation of \code{eigen(K)}).
#' 
#' The function \code{lmm.profile.restricted.likelihood} computes a profile restricted
#' likelihood: the values of \eqn{\tau}{tau} and \eqn{\sigma^2}{sigma^2} which
#' maximizes the likelihood are computed under the constraint
#' \eqn{ {\tau \over \tau + \sigma^2 } = h^2 }{tau/(tau + sigma^2) = h^2},
#' and the profiled likelihood value for these parameters is computed.
#' Whenever \eqn{k = 1}, it is similar to \code{lmm.diago.likelihood(h2 = h2, Y = Y, X = X, eigenK = eigen(K))}.
#' 
#' 
#' 
#' 
#' @return
#' The restricted likelihood value.
#' @seealso  \code{\link{lmm.diago.likelihood}}, \code{\link{lmm.diago}}, \code{\link{lmm.aireml}}
#' 
#' 
#' @keywords  Eigen decomposition   Likelihood
#' @examples
#' 
#' # Load data
#' data(AGT)
#' x <- as.bed.matrix(AGT.gen, AGT.fam, AGT.bim)
#' 
#' # Compute Genetic Relationship Matrix and its eigen decomposition
#' K <- GRM(x)
#' eiK <- eigen(K)
#' 
#' # simulate a phenotype
#' set.seed(1)
#' y <- 1 + lmm.simu(tau = 1, sigma2 = 2, eigenK = eiK)$y
#' 
#' # compute restricted likelihood for tau = 0.2 and s2 = 0.8
#' lmm.restricted.likelihood(y, K=K, tau = 0.2, s2 = 0.8)
#' 
#' # compute profile restricted likelihood for h2 = 0.2
#' lmm.profile.restricted.likelihood(y, K=K, h2 = 0.2)
#' 
#' # identity with the values computed with the diagonalisation trick
#' lmm.diago.likelihood(tau = 0.2, s2 = 0.8, Y = y, eigenK = eiK)
#' lmm.diago.likelihood(h2 = 0.2, Y = y, eigenK = eiK)
#' 
#' @export lmm.restricted.likelihood
lmm.restricted.likelihood <- function(Y, X = matrix(1, nrow = length(Y)), K, tau, s2) {
  if(!is.matrix(X) & !is.null(X)) stop("X should be a matrix or NULL")
  if(!is.vector(Y) & !is.matrix(Y)) 
    stop("Y should be a vector or a one-column matrix");
  if(is.matrix(Y)) {
    if(ncol(Y)!=1) 
      stop("Y should be a vector or a one-column matrix");
  } 
  if(is.matrix(K)) K <- list(K)
  
  if( any(is.na(Y)) ) {
    w <- !is.na(Y)
    if (!is.null(X)) X <- as.matrix(X[w,])
    Y <- Y[w]
    K <- lapply(K, function(x) x[w,w])
    warning(sum(!w), 'missing values are ignored.\n')
  }

  n <- length(Y)
  theta <- c(s2, tau)

  # X = NULL pour supprimer les effets fixes, y compris l'intercept
  if(is.null(X)) {
    if(any(sapply(K,nrow) != n) | any(sapply(K,ncol) != n))
      stop("Dimensions of Y and K mismatch")
    return( .Call(`_gaston_re_likelihood_nofix`, PACKAGE = "gaston", Y, K, theta) );
  }

  # sinon, X = matrice d'effets fixes
  if(nrow(X) != n) stop("Dimensions of X and Y mismatch")
  if(ncol(X) >= n) stop("Too many columns in X")
  if(any(sapply(K,nrow) != n) | any(sapply(K,ncol) != n))
    stop("Dimensions of Y and K mismatch")
  return( .Call(`_gaston_re_likelihood`, PACKAGE = "gaston", Y, X, K, theta) );
}

