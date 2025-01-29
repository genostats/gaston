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

