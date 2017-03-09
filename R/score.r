### P value with Davies method
davies <- function(q,lambda,h = rep(1,length(lambda)),delta = rep(0,length(lambda)),sigma=0,lim=10000,acc=0.0001) {
  r <- length(lambda)
  if (length(h) != r) stop("lambda and h should have the same length!")
  if (length(delta) != r) stop("lambda and delta should have the same length!")
  
  out <- .C("qfc",lambdas=as.double(lambda),noncentral=as.double(delta),df=as.integer(h),r=as.integer(r),sigma=as.double(sigma),q=as.double(q),lim=as.integer(lim),
             acc=as.double(acc),trace=as.double(rep(0,7)),ifault=as.integer(0),res=as.double(0), PACKAGE = "gaston")
  
  if (out$ifault==1) warning('In Davies method : Requested accuracy could not be obtained.')
  if (out$ifault==2) warning('In Davies method : Round-off error possibly significant.')
  if (out$ifault==3) stop('In Davies method : Invalid parameters.')
  if (out$ifault==4) stop('In Davies method : Unable to locate integration parameters.')
  
  out$res <- 1 - out$res
  if (out$res<0) out$res <- 0 
  return(out$res)
  
}


#### Test of variance parameter
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

score.variance.logistic <- function(K0, Y, X = matrix(1, length(Y)), K=NULL, acc_davies=1e-10, ...) {
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
    if(length(Y) != nrow(X)) stop("Dimensions of Y and X mismatch")
    
    if ( any( rowSums(is.na(X))>0 ) )
    {
      w <- ( rowSums(is.na(X))==0 )
      X <- as.matrix(X[w,])
      Y <- Y[w]
      if (is.list(K)) K <- lapply(K, function(x) x[w,w]) else if (is.matrix(K)) K <- K[w,w]
      K0 <- K0[w,w]
      warning(sum(!w), ' individuals with missing covariates are ignored.\n')
    } }

  if (is.list(K) | is.matrix(K))
  {
    model <- logistic.mm.aireml(Y, X=X, K=K, get.P = TRUE, ...)
    omega <- model$BLUP_omega
    if (!is.null(X)) omega <- omega + X%*%model$BLUP_beta
    pi <- 1/(1+exp(-omega))	
 
    eig <- eigen(model$P)
    eig$values[ eig$values<0 ] <- 0
    PP <- eig$vectors%*%( sqrt(eig$values)*t(eig$vectors) )
    T <- PP%*%K0%*%PP/2
    s <- t(Y-pi)%*%K0%*%(Y-pi)/2
  } else {
    if (!is.null(X))
    {
      model <- glm(Y~X-1, family=binomial())
      pi <- model$fitted.values
      V <- diag( sqrt(pi*(1-pi)) )
      P <- V - V%*%X%*%solve(t(X)%*%(V**2)%*%X)%*%t(X)%*%V**2
    } else {
      pi <- rep(1/2, length(Y)) 
      V <- diag( sqrt(pi*(1-pi)) )
      P <- V
    }
    T <- P%*%K0%*%t(P)/2
    s <- t(Y-pi)%*%K0%*%(Y-pi)/2
  }
  l <- eigen(T, symmetric=TRUE, only.values=TRUE)$values
  l <- l[l>0]
  #l <- l[l>mean(l)*1e-5]
  p <- davies(s, lambda = l, acc=acc_davies)
  
  return(list(score=s,p=p))
}


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


score.fixed.logistic <- function(x, Y, X = matrix(1,  length(Y)), K, ...) {
  if ( !is.matrix(K) & !is.list(K) ) stop("K must be a matrix, a list of matrix") 
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
  
  model <- logistic.mm.aireml(Y, X=X, K=K, get.P = TRUE, ...)
    
  omega <- model$BLUP_omega
  if (!is.null(X)) omega <- omega + X%*%model$BLUP_beta
  pi <- 1/(1+exp(-omega))	
  
  T <- t(x)%*%(Y-pi)
  V <- t(x)%*%model$P%*%x
  
  s <- t(T)%*%solve(V)%*%T
  logp <- pchisq( s, df = ncol(x), lower.tail=FALSE, log.p=TRUE)
  
  return(list(score=s,p=exp(logp), log.p=logp))
}
