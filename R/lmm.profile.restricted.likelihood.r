lmm.profile.restricted.likelihood <- function(Y, X = matrix(1, nrow = length(Y)), K, h2) {
  if(any(h2 < 0) | sum(h2) > 1) return(NA)

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

  # X = NULL pour supprimer les effets fixes, y compris l'intercept
  if(is.null(X)) {
    if(any(sapply(K,nrow) != n) | any(sapply(K,ncol) != n))
      stop("Dimensions of Y and K mismatch")
    return( .Call("gg_pre_likelihood_nofix", PACKAGE = "gaston", Y, K, h2) );
  }

  # sinon, X = matrice d'effets fixes
  if(nrow(X) != n) stop("Dimensions of X and Y mismatch")
  if(ncol(X) >= n) stop("Too many columns in X")
  if(any(sapply(K,nrow) != n) | any(sapply(K,ncol) != n))
    stop("Dimensions of Y and K mismatch")
  return( .Call("gg_pre_likelihood", PACKAGE = "gaston", Y, X, K, h2) );
}

