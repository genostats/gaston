#' Linear mixed model data simulation
#' 
#' @description
#' Simulate data under a linear mixed model, using the eigen decomposition
#' of the variance matrix.
#' 
#' @param tau  Model parameter
#' @param sigma2  Model parameter
#' @param K  (Optional) A positive symmetric matrix \eqn{K}
#' @param eigenK  Eigen decomposition of \eqn{K}
#' @param X  Covariable matrix
#' @param beta  Fixed effect vector of covariables
#' 
#' @details
#' The data are simulated under the following linear mixed model :
#' \deqn{ Y = X\beta + \omega + \varepsilon }{ Y = X beta + omega + epsilon }
#' with \eqn{ \omega \sim N(0,\tau K) }{omega ~ N(0, tau K)} and
#' \eqn{ \varepsilon \sim N(0,\sigma^2 I_n) }{epsilon ~ N(0, sigma^2 I_n)}.
#' 
#' The simulation uses \eqn{K} only through its eigen decomposition; the parameter
#' \code{K} is therefore optional.
#' 
#' @return
#' A named list with two members:
#' \item{y}{ Simulated value of \eqn{Y} }
#' \item{omega}{Simulated value of \eqn{\omega}{omega}  }
#' @seealso  \code{\link{random.pm}}
#' 
#' @keywords  Simulations
#' @examples
#' 
#' # generate a random positive matrix 
#' set.seed(1)
#' R <- random.pm(503)
#' 
#' # simulate data with a "polygenic component" 
#' y <-  lmm.simu(0.3, 1, eigenK = R$eigen)
#' str(y)
#' 
#' @export lmm.simu
lmm.simu <- function(tau, sigma2, K, eigenK = eigen(K), X, beta) {
  if(any(eigenK$values < -1e-5))
    warning("K is not positive, setting negative eigenvalues to 0")
  eigenK$values[eigenK$values < 0] <- 0
  # partie fixe
  xbeta <- if(!missing(X)) {
    if(nrow(X) != nrow(eigenK$vectors)) stop("Dimensions mismatch (K and X)")
    if(length(beta) != ncol(X)) stop("Dimensions mismatch (X and beta)")
    X %*% beta
  } else 0
  # partie aléatoire
  omega <- eigenK$vectors %*% rnorm( nrow(eigenK$vectors), sd = sqrt(tau*eigenK$values) ) 
  y <- xbeta + omega + rnorm( nrow(eigenK$vectors), sd=sqrt(sigma2) )
  return( list(y = y, omega = omega) )
}
