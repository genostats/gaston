#' Logistic mixed model fitting with Penalized Quasi-Likelihood / AIREML
#' 
#' @description
#' Estimate the parameters of a logistic linear mixed model using the Penalized
#' Quasi-Likelihood with an AIREML step for the linear model.
#' 
#' @param Y  Binary phenotype vector
#' @param X  Covariable matrix. By default, a column of ones to include an intercept in the model
#' @param K  A positive definite matrix or a \code{list} of such matrices
#' @param min_tau  Minimal value for model parameter \eqn{\tau}{tau} (if missing, will be set to \eqn{10^{-6}}{1e-6})
#' @param tau  (Optional) Optimization starting point for variance component(s) \code{tau}
#' @param beta  (Optional) Optimization starting point for fixed effect(s) \code{beta}
#' @param constraint  If \code{TRUE}, the model parameters respect the contraints given by \code{min_tau}
#' @param max.iter  Maximum number of iterations
#' @param eps  The algorithm stops when the gradient norm is lower than this parameter
#' @param verbose  If \code{TRUE}, display information on the algorithm progress
#' @param get.P If \code{TRUE}, the function sends back the last matrix \eqn{P} computed in the optimization process
#' @param EM If \code{TRUE}, the AIREML step is replaced by an EM step
#' 
#' 
#' @details
#' Estimate the parameters of the following logistic mixed model:
#' \deqn{ \mbox{logit}(P[Y=1|X,\omega_1,\ldots,\omega_k])  = X\beta + \omega_1 + \ldots + \omega_k}{logit P(Y=1|X,omega_1,...,omega_k)  = X beta  + omega_1 + ... + omega_k}
#' with \eqn{ \omega_i \sim N(0,\tau_i K_i) }{omega_i ~ N(0, tau_i K_i)} for \eqn{ i \in 1, \dots,k }.
#' 
#' The estimation is based on the Penalized Quasi-Likelihood with an AIREML step for the linear model
#' (the algorithm is similar to the algorithm described in Chen et al 2016). If \code{EM = TRUE}
#' the AIREML step is replaced by an EM step. In this case the convergence will be much slower,
#' you're advised to use a large value of \code{max.iter}.
#' 
#' The variance matrices \eqn{K_1}{K_1}, ..., \eqn{K_k}{K_k}, are specified through the parameter \code{K}.
#' 
#' After convergence, the function also compute Best Linear Unbiased Predictors (BLUPs) for
#' \eqn{\beta}{beta} and \eqn{\omega}{omega}.
#' 
#' 
#' 
#' @return
#' A named list with members:
#' \item{tau}{ Estimate(s) of the model parameter(s) \eqn{\tau_1, \dots, \tau_k}{tau_1, ..., tau_k} }
#' \item{niter}{ Number of iterations done }
#' \item{P}{ Last computed value of matrix P (see reference) }
#' \item{BLUP_omega}{ BLUPs of random effects }
#' \item{BLUP_beta}{ BLUPs of fixed effects \eqn{\beta}{beta}}
#' \item{varbeta}{ Variance matrix for \eqn{\beta}{beta} estimates }
#' If \code{get.P = TRUE}, there is an additional member:
#' \item{P}{The last matrix \eqn{P} computed in the AIREML step}
#' @seealso  \code{\link{lmm.aireml}}, \code{\link{lmm.diago}}, \code{\link{lmm.simu}}
#' 
#' @references Gilmour, A. R., Thompson, R., & Cullis, B. R. (1995), \emph{Average information REML: an efficient algorithm for variance
#' parameter estimation in linear mixed models}, Biometrics, \bold{1440-1450}
#' 
#' Chen, Han et al. (2016), \emph{Control for Population Structure and Relatedness for Binary Traits in Genetic Association Studies via
#' Logistic Mixed Models}, The American Journal of Human Genetics, \bold{653--666}
#' 
#' @keywords  Penalized Quasi-Likelihood   Average Information Restricted Maximum Likelihood (AIREML)   Genetic effect
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
#' mu <- 1 + x %*% rnorm(ncol(x), sd = 2)/sqrt(ncol(x))
#' pi <- 1/(1+exp(-mu))
#' y <- 1*( runif(length(pi))<pi )
#' 
#' # Estimates
#' estimates <- logistic.mm.aireml(y, K = K, verbose = FALSE)
#' str(estimates)
#' 
#' @export logistic.mm.aireml
logistic.mm.aireml <- function(Y, X = matrix(1, nrow = length(Y)), K, min_tau, tau, beta, constraint = TRUE, max.iter = 50L, eps = 1e-5, 
                           verbose = getOption("gaston.verbose",TRUE), get.P = FALSE, EM = FALSE) {

  if(!is.vector(Y) & !is.matrix(Y)) stop("Y should be a vector or a one-column matrix")
  if(is.matrix(Y)) { if(ncol(Y)!=1) stop("Y should be a vector or a one-column matrix")  } 
  if(!is.matrix(X) & !is.null(X)) stop("X should be a matrix or NULL")
  
  if( any(is.na(Y)) ) {
    w <- !is.na(Y)
    if (!is.null(X)) X <- as.matrix(X[w,])
    Y <- Y[w]
    if (is.list(K)) K <- lapply(K, function(x) x[w,w])
    if (is.matrix(K)) K <- K[w,w]
    warning(sum(!w), "missing phenotype values are ignored.\n")
  }
  
  if (sum(Y %in% c(0,1))<length(Y)) stop("Y should contain only '0' and '1' values")   
  n <- length(Y)
 
   
  # on s'occupe de tau et start_tau
  if(is.matrix(K)) {
    if(missing(tau)) {
      tau <- 0.1;
      start_tau <- TRUE      # !!! modifié par rv juillet 2018 : le tau calculé dans le code C++ ne marche pas bien
    } else start_tau <- TRUE
  } else if(is.list(K)) {
    if(missing(tau)) {
      tau <- rep(0.1, length(K))
      start_tau <- TRUE      # !!! idem ci-dessus
    } else start_tau <- TRUE
  } else stop("K should be a matrix or a list of matrices");
  
  # X = NULL pour supprimer les effets fixes, y compris l'intercept
  if(is.null(X)) {
    if(!missing(beta)) warning("No covariable matrix X, parameter 'beta' will be ignored \n")
    if(is.matrix(K)) {
      if(nrow(K) != n | ncol(K) != n) stop("Dimensions of Y and K mismatch")
      if(missing(min_tau)) min_tau <- 1e-6
      return( .Call(`_gaston_AIREML1_logit_nofix`, PACKAGE = "gaston", Y, K, constraint, min_tau, max.iter, eps, verbose, tau, start_tau, get.P, EM) )
    } 
    else if(is.list(K)) {
      if(any(sapply(K,nrow) != n) | any(sapply(K,ncol) != n))
        stop("Dimensions of Y and K mismatch")
      if(missing(min_tau)) min_tau <- rep(1e-6, length(K))
      return( .Call(`_gaston_AIREMLn_logit_nofix`, PACKAGE = "gaston", Y, K, constraint, min_tau, max.iter, eps, verbose, tau, start_tau, get.P, EM) )
    }
  }

  # sinon, X = matrice d'effets fixes
  if(nrow(X) != n) stop("Dimensions of X and Y mismatch")
  
  if( max( rowSums(is.na(X)) )>0 ) {
    w <- rowSums(is.na(X))==0
    X <- as.matrix(X[w,])
    Y <- Y[w]
    if (is.list(K)) K <- lapply(K, function(x) x[w,w]) else K <- K[w,w]
    warning(sum(!w), " individuals with missing covariates are ignored.\n")
  }
  
  n <- length(Y)
  if(ncol(X) >= n) stop("Too many columns in X")

  # if(missing(beta)) {
  #  beta <- glm(Y~X-1, family=binomial)$coefficients
  # } else { if (length(beta)!=ncol(X)) stop("Dimensions of X and beta mismatch") }
  # start_beta <- TRUE
  if(missing(beta)) {
    beta <- rep(0, ncol(X))
    start_beta <- FALSE
  } else { 
    if(length(beta) != ncol(X)) 
      stop("Dimensions of X and beta mismatch")
    start_beta <- TRUE
  }

  if(is.matrix(K)) {
    if(nrow(K) != n | ncol(K) != n) stop("Dimensions of Y and K mismatch")
    if(missing(min_tau)) min_tau <- 1e-6
    return( .Call(`_gaston_AIREML1_logit`, PACKAGE = "gaston", Y, X, K, constraint, min_tau, max.iter, eps, verbose, tau, beta, start_tau, start_beta, get.P, EM) )
  }
  else if(is.list(K)) {
    if(any(sapply(K,nrow) != n) | any(sapply(K,ncol) != n))
      stop("Dimensions of Y and K mismatch")
    if(missing(min_tau)) min_tau <- rep(1e-6, length(K))
    return( .Call(`_gaston_AIREMLn_logit`, PACKAGE = "gaston", Y, X, K, constraint, min_tau, max.iter, eps, verbose, tau, beta, start_tau, start_beta, get.P, EM) )
  }
}

