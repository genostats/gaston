#' @rdname score.fixed.logistic
#' @export
score.fixed.linear <- function(x, Y, X = matrix(1, length(Y)), K, ...) {
  if ( !is.matrix(K) & !is.list(K) ) stop("K must be a matrix or a list of matrix") 
  if (!is.matrix(x)) stop("'x' must to be a matrix") 
  if(length(Y) != nrow(x)) stop("Dimensions of Y and x mismatch")
  if (is.matrix(K)) if ( length(Y) != nrow(K) | length(Y) != ncol(K) ) stop("Dimensions of Y and K mismatch")
  if (is.list(K)) if (length( unique(unlist(lapply(K, function(x) c(nrow(x),ncol(x)) )) ))!=1 ) stop("Matrix in K must be square and of same dimensions")
  if (is.list(K)) if  (length(Y) != unique(unlist(lapply(K, function(x) c(nrow(x),ncol(x)) )) )) stop("Dimensions of Y and K mismatch")
  if(!is.matrix(X) & !is.null(X)) stop("X should be a matrix or NULL")
  
  if( any(is.na(Y)) | any( rowSums(is.na(x))>0 ) )
  {
    w <- ( !is.na(Y) & rowSums(is.na(x))==0 )
    if (!is.null(X)) X <- as.matrix(X[w,])
    Y <- Y[w]
    if (is.list(K)) K <- lapply(K, function(x) x[w,w]) else K <- K[w,w]
    x <- as.matrix(x[w,])
    warning(sum(!w), ' individuals with missing phenotype are ignored.\n')
  }
  
  if (!is.null(X))
  {
    if(length(Y) != nrow(X)) stop("Dimensions of Y and X mismatch")
    
    if ( any( rowSums(is.na(X))>0 ) )
    {
      w <- ( rowSums(is.na(X))==0 )
      X <- as.matrix(X[w,])
      Y <- Y[w]
      if (is.list(K)) K <- lapply(K, function(x) x[w,w]) else if (is.matrix(K)) K <- K[w,w]
      x <- as.matrix(x[w,])
      warning(sum(!w), ' individuals with missing covariates are ignored.\n')
   } }

  model <- lmm.aireml(Y, X = X, K = K, get.P = TRUE, ...)
  XtP <- t(x)%*%model$P
  V <- XtP%*%x 
  T <- XtP%*%Y 
  
  s <- t(T)%*%solve(V)%*%T
  logp <- pchisq( s, df = ncol(x), lower.tail=FALSE, log.p=TRUE)
  
  return(list(score=s,p=exp(logp), log.p=logp))
}
