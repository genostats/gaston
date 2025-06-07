#' Linear mixed model fitting with AIREML
#' 
#' @description
#' Estimate the parameters of a linear mixed model,
#' using Average Information Restricted Maximum Likelihood (AIREML) algorithm.
#' 
#' @param Y  Phenotype vector
#' @param X  Covariable matrix. By default, a column of ones to include an intercept in the model
#' @param K  A positive definite matrix or a \code{list} of such matrices
#' @param EMsteps  Number of EM steps ran prior the AIREML
#' @param EMsteps_fail  Number of EM steps performed when the AIREML algorithm fail to improve the likelihood value
#' @param EM_alpha  Tweaking parameter for the EM (see Details)
#' @param min_tau  Minimal value for model parameter \eqn{\tau}{tau} (if missing, will be set to \eqn{10^{-6}}{1e-6})
#' @param min_s2  Minimal value for model parameter \eqn{\sigma^2}{sigma^2}
#' @param theta  (Optional) Optimization starting point \code{theta = c(sigma^2, tau)}
#' @param constraint  If \code{TRUE}, the model parameters respect the contraints given by \code{min_tau} and \code{min_s2}
#' @param max_iter  Maximum number of iterations
#' @param eps  The algorithm stops when the gradient norm is lower than this parameter
#' @param verbose  If \code{TRUE}, display information on the algorithm progress
#' @param contrast  If \code{TRUE}, use a contrast matrix to compute the Restricted Likelihood (usually slower)
#' @param get.P If \code{TRUE}, the function sends back the last matrix \eqn{P} computed in the optimization process
#' 
#' 
#' @details
#' Estimate the parameters of the following linear mixed model, using AIREML algorithm:
#' \deqn{ Y = X\beta + \omega_1 + \ldots + \omega_k + \varepsilon }{ Y =X beta + omega_1 + ... + omega_k + epsilon }
#' with \eqn{ \omega_i \sim N(0,\tau_i K_i) }{omega_i ~ N(0, tau_i K_i)} for \eqn{ i \in 1, \dots,k } and
#' \eqn{ \varepsilon \sim N(0,\sigma^2 I_n) }{epsilon ~ N(0, sigma^2 I_n)}.
#' 
#' The variance matrices \eqn{K_1}{K_1}, ..., \eqn{K_k}{K_k}, are specified through the parameter \code{K}.
#' 
#' If \code{EMsteps} is positive, the function will use this number of EM steps to compute a better starting point
#' for the AIREML algorithm. Setting \code{EMsteps} to a value higher than \code{max_iter} leads to an EM optimization.
#' It can happen that after an AIREML step, the likelihood did not increase: if this
#' happens, the functions falls back to \code{EMsteps_fail} EM steps. The parameter \code{EM_alpha} can be set to
#' a value higher than \eqn{1} to attempt to accelerate EM convergence; this could also result in uncontrolled
#' behaviour and should be used with care.
#' 
#' After convergence, the function also compute Best Linear Unbiased Predictors (BLUPs) for
#' \eqn{\beta}{beta} and \eqn{\omega}{omega}, and an
#' estimation of the participation of the fixed effects to the variance of \eqn{Y}.
#' 
#' @return
#' A named list with members:
#' \item{sigma2}{ Estimate of the model parameter \eqn{\sigma^2}{sigma^2} }
#' \item{tau}{ Estimate(s) of the model parameter(s) \eqn{\tau_1, \dots, \tau_k}{tau_1, ..., tau_k} }
#' \item{logL}{ Value of log-likelihood }
#' \item{logL0}{ Value of log-likelihood under the null model (without random effect) }
#' \item{niter}{ Number of iterations done }
#' \item{norm_grad}{ Last computed gradient's norm }
#' \item{P}{ Last computed value of matrix P (see reference) }
#' \item{Py}{ Last computed value of vector Py (see reference) }
#' \item{BLUP_omega}{ BLUPs of random effects }
#' \item{BLUP_beta}{ BLUPs of fixed effects \eqn{\beta}{beta}}
#' \item{varbeta}{ Variance matrix for \eqn{\beta}{beta} estimates }
#' \item{varXbeta}{ Participation of fixed effects to variance of Y }
#' If \code{get.P = TRUE}, there is an additional member:
#' \item{P}{The last matrix \eqn{P} computed in the AIREML step}
#' @seealso \code{\link{lmm.diago}}, \code{\link{logistic.mm.aireml}}, \code{\link{lmm.simu}}
#' 
#' @references Gilmour, A. R., Thompson, R., & Cullis, B. R. (1995), \emph{Average information REML: an efficient algorithm for variance parameter estimation in linear mixed models}, Biometrics, \bold{1440-1450}
#' 
#' @keywords  Average Information Restricted Maximum Likelihood (AIREML)   Genetic effect
#' @examples
#' 
#' # Load data
#' data(AGT)
#' x <- as.bed.matrix(AGT.gen, AGT.fam, AGT.bim)
#' 
#' # Compute Genetic Relationship Matrix
#' standardize(x) <- "p"
#' K <- GRM(x)
#' 
#' # Simulate a quantitative genotype under the LMM
#' set.seed(1)
#' y <- 1 + x %*% rnorm(ncol(x), sd = 1)/sqrt(ncol(x)) + rnorm(nrow(x), sd = sqrt(2))
#' 
#' # Estimates
#' estimates <- lmm.aireml(y, K = K, verbose = FALSE)
#' str(estimates)
#' 
#' @export lmm.aireml
lmm.aireml <- function(Y, X = matrix(1, nrow = length(Y)), K, EMsteps = 0L, EMsteps_fail = 1L, EM_alpha = 1, 
                   min_tau, min_s2 = 1e-6, theta, constraint = TRUE, max_iter = 50L, eps = 1e-5, 
                   verbose = getOption("gaston.verbose",TRUE), contrast = FALSE, get.P = FALSE) {
  if(!is.matrix(X) & !is.null(X)) stop("X should be a matrix or NULL")
  if(!is.vector(Y) & !is.matrix(Y)) 
    stop("Y should be a vector or a one-column matrix");
  if(is.matrix(Y)) {
    if(ncol(Y)!=1) 
      stop("Y should be a vector or a one-column matrix");
  } 
  
  if( any(is.na(Y)) ) {
    w <- !is.na(Y)
    if (!is.null(X)) X <- as.matrix(X[w,])
    Y <- Y[w]
    if (is.list(K)) K <- lapply(K, function(x) x[w,w]) else K <- K[w,w]
    warning(sum(!w), 'missing values are ignored.\n')
  }

  n <- length(Y)

  # on s'occupe de theta et start_theta
  if(is.matrix(K)) {
    if(missing(theta)) {
      theta <- c(0,0);
      start_theta <- FALSE
    } else start_theta <- TRUE
  } else if(is.list(K)) {
    if(missing(theta)) {
      theta <- rep(0, length(K)+1)
      start_theta <- FALSE
    } else start_theta <- TRUE
  } else stop("K should be a matrix or a list of matrices");

  # X = NULL pour supprimer les effets fixes, y compris l'intercept
  if(is.null(X)) {
    if(is.matrix(K)) {
      if(nrow(K) != n | ncol(K) != n) stop("Dimensions of Y and K mismatch")
      if(missing(min_tau)) min_tau <- 1e-6
      return( .Call(`_gaston_AIREML1_nofix`, PACKAGE = "gaston", Y, K, EMsteps, EMsteps_fail, EM_alpha, constraint, 
                    min_s2, min_tau, max_iter, eps, verbose, theta, start_theta, get.P) )
    } 
    else if(is.list(K)) {
      if(any(sapply(K,nrow) != n) | any(sapply(K,ncol) != n))
        stop("Dimensions of Y and K mismatch")
      if(missing(min_tau)) min_tau <- rep(1e-6, length(K))
      return( .Call(`_gaston_AIREMLn_nofix`, PACKAGE = "gaston", Y, K, EMsteps, EMsteps_fail, EM_alpha, constraint, 
                    min_s2, min_tau, max_iter, eps, verbose, theta, start_theta, get.P) )
    }
  }

  # sinon, X = matrice d'effets fixes
  if(nrow(X) != n) stop("Dimensions of X and Y mismatch")
  if(ncol(X) >= n) stop("Too many columns in X")
  if(is.matrix(K)) { # only on K
    if(nrow(K) != n | ncol(K) != n) stop("Dimensions of Y and K mismatch")
    if(missing(min_tau)) min_tau <- 1e-6
    if(contrast) { 
      return( .Call(`_gaston_AIREML1_contrast`, PACKAGE = "gaston", Y, X, K, EMsteps, EMsteps_fail, EM_alpha, constraint, min_s2, min_tau, 
                  max_iter, eps, verbose, theta, start_theta, get.P) )
    } else {
      return( .Call(`_gaston_AIREML1`, PACKAGE = "gaston", Y, X, K, EMsteps, EMsteps_fail, EM_alpha, constraint, min_s2, min_tau, 
                  max_iter, eps, verbose, theta, start_theta, get.P) )
    }
  } else if(is.list(K)) { # many K's
    if(any(sapply(K,nrow) != n) | any(sapply(K,ncol) != n))
      stop("Dimensions of Y and K mismatch")
    if(missing(min_tau)) min_tau <- rep(1e-6, length(K))
    if(contrast) {
      return( .Call(`_gaston_AIREMLn_contrast`, PACKAGE = "gaston", Y, X, K, EMsteps, EMsteps_fail, EM_alpha, constraint, min_s2, min_tau, 
                  max_iter, eps, verbose, theta, start_theta, get.P) )
    } else {
      return( .Call(`_gaston_AIREMLn`, PACKAGE = "gaston", Y, X, K, EMsteps, EMsteps_fail, EM_alpha, constraint, min_s2, min_tau, 
                  max_iter, eps, verbose, theta, start_theta, get.P) )
    }
  }
}

