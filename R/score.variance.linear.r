#' @rdname score.variance.logistic
#' @export
score.variance.linear <- function(K0, Y, X = matrix(1, length(Y)), K=NULL, acc_davies=1e-10, ...) {
  if ( !is.null(K) & !is.matrix(K) & !is.list(K) ) stop("K must be a matrix, a list of matrix or 'NULL'") 
  if (!is.matrix(K0)) stop("'K0' must to be a matrix") 
  if(length(Y) != nrow(K0) | length(Y) != ncol(K0)  ) stop("Dimensions of Y and K0 mismatch")
  if (is.matrix(K)) if ( length(Y) != nrow(K) | length(Y) != ncol(K) ) stop("Dimensions of Y and K mismatch")
  if (is.list(K)) if (length( unique(unlist(lapply(K, function(x) c(nrow(x),ncol(x)) )) ))!=1 ) stop("Matrix in K must be square and of same dimensions")
  if (is.list(K)) if  (length(Y) != unique(unlist(lapply(K, function(x) c(nrow(x),ncol(x)) )) )) stop("Dimensions of Y and K mismatch")
  if(!is.matrix(X) & !is.null(X)) stop("X should be a matrix or NULL")
  
  if( any(is.na(Y)) )
  {
    w <- !is.na(Y)
    if (!is.null(X)) X <- as.matrix(X[w,])
    Y <- Y[w]
    if (is.list(K)) K <- lapply(K, function(x) x[w,w]) else if (is.matrix(K)) K <- K[w,w]
    K0 <- K0[w,w]
    warning(sum(!w), ' individuals with missing phenotype are ignored.\n')
  }
  
  if (!is.null(X))
  {
    r <- ncol(X)
    if(length(Y) != nrow(X)) stop("Dimensions of Y and X mismatch")
    
    if ( any( rowSums(is.na(X))>0 ) )
    {
      w <- (rowSums(is.na(X))==0)
      X <- as.matrix(X[w,])
      Y <- Y[w]
      if (is.list(K)) K <- lapply(K, function(x) x[w,w]) else if (is.matrix(K)) K <- K[w,w]
      K0 <- K0[w,w]
      warning(sum(!w), ' individuals with missing covariates are ignored.\n')
    }
  } else r <- 0 
    
  if (is.list(K) | is.matrix(K))
  {
    model <- lmm.aireml(Y, X = X, K = K, get.P = TRUE, ...)
    eig <- eigen(model$P)
    eig$values[ eig$values<0 ] <- 0
    PP <- eig$vectors%*%( sqrt(eig$values)*t(eig$vectors) )
    T <- PP%*%K0%*%PP/2
    s <- (t(Y)%*%model$P)%*%K0%*%(model$P%*%Y)/2
  } else {
    if (!is.null(X)) P <- diag(1, length(Y)) - X %*% solve(t(X) %*% X) %*%t(X) else P <- diag(1, length(Y))
    sigma2 <- (sum((P %*% Y)^2)/(length(Y) - r))
    T <- P %*% K0 %*% P /2
    s <- t(Y) %*% T %*% Y
    T <- sigma2 * T
  }  
  l <- eigen(T, symmetric=TRUE, only.values=TRUE)$values
  l <- l[l>0]
  #l <- l[l>mean(l)*1e-5]
  p <- davies(s, lambda = l, acc=acc_davies)
  
  return(list(score=s,p=p))
}

