#' Score Test for Covariates with Fixed Effects in Linear or Logistic Mixed Model
#'
#' @description  Score Test for association between covariates and phenotype.
#' 
#' @param x  A matrix of covariates
#' @param Y  The phenotype vector
#' @param X  A covariable matrix. The default is a column vector of ones, to include an intercept in the model
#' @param K  A positive definite matrix or a \code{list} of such matrices
#' @param ...  Optional arguments used to fit null model in \code{lmm.aireml} or \code{logistic.mm.aireml} function.
#' 
#' 
#' @details
#' The function \code{score.fixed.linear} considers the linear mixed model
#' \deqn{ Y = X\alpha + x\beta + \omega_1 + \ldots + \omega_k + \varepsilon }{ Y = X alpha + x beta + omega_1 + ... + omega_k + varepsilon }
#' whereas the \code{score.fixed.logistic} function considers the following logistic model
#' \deqn{ \mbox{logit}(P[Y=1|X,x,\omega_1,\ldots,\omega_k])  = X\alpha + x\beta + \omega_1 + \ldots + \omega_k}{logit(P[Y=1|X,x,omega_1,...,omega_k])  = Xalpha + x beta + omega_1 + ... + omega_k}
#' with \eqn{ \omega_j \sim N(0,\tau_j K_j) }{omega_j ~ N(0, tau_j K_j)} where \eqn{K_j} are Genetic Relationship Matrix (GRM), \eqn{ \varepsilon \sim N(0,\sigma^2 I_n) }{epsilon ~ N(0, sigma^2 I_n)}
#' and fixed effects \eqn{\alpha}{alpha} and \eqn{\beta}{beta}.
#' 
#' The two functions give score test for
#' \eqn{H_0}{H0} : \eqn{\beta=0}{beta=0} vs \eqn{H_1}{H1} : \eqn{\beta\neq 0}{beta !=0}.
#' In this aim, all parameters under null model are estimated with \code{lmm.aireml} or \code{logistic.mm.aireml}.
#' 
#' @return
#' A named list of values:
#' \item{score}{ Estimated score }
#' \item{p}{ The corresponding p-value }
#' \item{log.p}{ The logarithm of corresponding p-value }
#' @seealso  \code{\link{lmm.aireml}}, \code{\link{logistic.mm.aireml}}
#' 
#' 
#' @keywords  Score Test
#' @examples
#' # Load data
#' data(AGT)
#' x <- as.bed.matrix(AGT.gen, AGT.fam, AGT.bim)
#' standardize(x) <- "p"
#' 
#' # Calculate GRM et its eigen decomposition
#' k <- GRM(x)
#' eig <- eigen(k)
#' eig$values <- round(eig$values, 5)
#' 
#' # generate covariate matrix
#' set.seed(1)
#' X <- cbind( rbinom(nrow(x), 1, prob=1/2), rnorm(nrow(x)) )
#' 
#' 
#' # simulate quantitative phenotype with polygenic component and covariate effects
#' y <- X %*% c(-1,0.5) + lmm.simu(0.3,1,eigenK=eig)$y
#' 
#' t <- score.fixed.linear(X, y, K=k, verbose=FALSE)
#' str(t)
#' 
#' 
#' # simulate binary phenotype with polygenic component and covariate effects
#' mu <- X %*% c(-1,0.5) + lmm.simu(1, 0, eigenK=eig)$y
#' pi <- 1/(1+exp(-mu))
#' y <- 1*( runif(length(pi))<pi )
#' 
#' tt <- score.fixed.logistic(X, y, K=k, verbose=FALSE)
#' str(tt)
#' 
#' @export 
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
