#' Linear mixed model fitting with the diagonalization trick
#' 
#' @description
#' Estimate the parameters of a linear mixed model, using the "diagonalization trick".
#' 
#' @param Y  Phenotype vector
#' @param X  Covariable matrix
#' @param eigenK  Eigen decomposition of \eqn{K} (a positive symmetric matrix)
#' @param p  Number of Principal Components included in the mixed model with fixed effect
#' @param method  Optimization method to use
#' @param min_h2  Minimum admissible value
#' @param max_h2  Maximum admissible value
#' @param verbose  If \code{TRUE}, display information on the function actions
#' @param tol  Accuracy of estimation
#' 
#' 
#' @details
#' Estimate the parameters of the following linear mixed model, computing the restricted likelihood as in \code{lmm.diago.likelihood},
#' and using either a Newton algorithm, or Brent algorithm as in \code{\link[stats:optimize]{optimize}}:
#' \deqn{ Y = (X|PC)\beta + \omega + \varepsilon }{ Y = (X|PC) beta + omega + epsilon }
#' with \eqn{ \omega \sim N(0,\tau K) }{omega ~ N(0, tau K)} and
#' \eqn{ \varepsilon \sim N(0,\sigma^2 I_n) }{epsilon ~ N(0, sigma^2 I_n)}.
#' 
#' The matrix \eqn{K} is given through its eigen decomposition, as produced by \code{eigenK = eigen(K, symmetric = TRUE)}.
#' The matrix \eqn{(X|PC)} is the concatenation of the covariable matrix \eqn{X} and
#' of the first \eqn{p} eigenvectors of \eqn{K}, included in the model with fixed effects.
#' 
#' 
#' 
#' @return
#' If the parameter \code{p} is a scalar, a list with following elements :
#' \item{sigma2}{ Estimate of the model parameter \eqn{\sigma^2}{sigma^2} }
#' \item{tau}{ Estimate(s) of the model parameter(s) \eqn{\tau_1, \dots, \tau_k}{tau_1, ..., tau_k} }
#' \item{Py}{ Last computed value of vector Py (see reference) }
#' \item{BLUP_omega}{ BLUPs of random effects }
#' \item{BLUP_beta}{ BLUPs of fixed effects \eqn{\beta}{beta} (only the components corresponding to \eqn{X})}
#' \item{Xbeta}{ Estimate of \eqn{(X|PC)\beta}{(X|PC)beta} }
#' \item{varbeta}{ Variance matrix for \eqn{\beta}{beta} estimates (only the components corresponding to \eqn{X}) }
#' \item{varXbeta}{ Participation of fixed effects to variance of Y }
#' \item{p}{ Number of Principal Components included in the linear mixed model with fixed effect }
#' 
#' If the paramer \code{p} is a vector of length \code{> 1}, a \code{list} of lists as described above,
#' one for each value in \code{p}.
#' @seealso  \code{\link{lmm.diago.likelihood}}, \code{\link{lmm.aireml}}, \code{\link[stats:optimize]{optimize}}
#' 
#' 
#' @keywords  Eigen decomposition   Likelihood Maximization
#' @examples
#' 
#' # Load data
#' data(AGT)
#' x <- as.bed.matrix(AGT.gen, AGT.fam, AGT.bim)
#' 
#' # Compute Genetic Relationship Matrix
#' K <- GRM(x)
#' 
#' # eigen decomposition of K
#' eiK <- eigen(K)
#' 
#' # simulate a phenotype
#' set.seed(1)
#' y <- 1 + lmm.simu(tau = 1, sigma2 = 2, eigenK = eiK)$y
#' 
#' # Estimations
#' R <- lmm.diago(Y = y, eigenK = eiK, p = c(0,10))
#' str(R)
#' 
#' @export lmm.diago
lmm.diago <- function(Y, X = matrix(1, nrow=length(Y)), eigenK, p = 0, method = c("newton", "brent"), min_h2 = 0, max_h2 = 1, 
                      verbose = getOption("gaston.verbose", TRUE), tol = .Machine$double.eps^0.25) {
  if( any(is.na(Y)) ) 
    stop('Missing data in Y.')

  if(!is.null(X)) {
    if( length(Y)<(ncol(X)+max(p)) ) 
      stop('The total number of covariables and PCs as fixed effects cannot exceed the number of observations.')
    if( length(Y)!=nrow(X) ) 
      stop('Length of Y and the number of rows of X differ.')
  } else if( length(Y) < max(p) )
    stop('The total number of covariables and PCs as fixed effects cannot exceed the number of observations.')

  Sigma <- eigenK$values
  if( length(Y)!=length(Sigma) ) 
    stop('Length of Y and number of eigenvalues differ.')

  w <- which(Sigma < 1e-6)
  Sigma[w] <- 1e-6

  U <- eigenK$vectors
 
  method <- match.arg(method)
  brent <- (method == "brent") 
  if(is.null(X))
    return(.Call(`_gaston_fit_diago_nocovar`, PACKAGE = "gaston", Y, p, Sigma, U, min_h2, max_h2, tol, verbose, brent))
  else
    return(.Call(`_gaston_fit_diago`, PACKAGE = "gaston", Y, X, p, Sigma, U, min_h2, max_h2, tol, verbose, brent))
  
}

