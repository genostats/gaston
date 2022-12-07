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
      return( .Call("gg_AIREML1_logit_nofix",  PACKAGE = "gaston", Y, K, constraint, min_tau, max.iter, eps, verbose, tau, start_tau, get.P, EM) )
    } 
    else if(is.list(K)) {
      if(any(sapply(K,nrow) != n) | any(sapply(K,ncol) != n))
        stop("Dimensions of Y and K mismatch")
      if(missing(min_tau)) min_tau <- rep(1e-6, length(K))
      return( .Call("gg_AIREMLn_logit_nofix", PACKAGE = "gaston", Y, K, constraint, min_tau, max.iter, eps, verbose, tau, start_tau, get.P, EM) )
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
    return( .Call("gg_AIREML1_logit",  PACKAGE = "gaston", Y, X, K, constraint, min_tau, max.iter, eps, verbose, tau, beta, start_tau, start_beta, get.P, EM) )
  }
  else if(is.list(K)) {
    if(any(sapply(K,nrow) != n) | any(sapply(K,ncol) != n))
      stop("Dimensions of Y and K mismatch")
    if(missing(min_tau)) min_tau <- rep(1e-6, length(K))
    return( .Call("gg_AIREMLn_logit", PACKAGE = "gaston", Y, X, K, constraint, min_tau, max.iter, eps, verbose, tau, beta, start_tau, start_beta, get.P, EM) )
  }
}

