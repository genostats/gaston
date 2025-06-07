#' Variance Component Test in Linear or Logistic Mixed Model
#' 
#' @description  Test if a variance component is significaly different from 0 using score test in a Linear or Logistic Mixed Model.
#' 
#' @param K0  A positive definite matrix
#' @param Y  The phenotype vector
#' @param X  A covariate matrix. The default is a column vector of ones, to include an intercept in the model
#' @param K  A positive definite matrix or a \code{list} of such matrices
#' @param acc_davies  Accuracy in Davies method used to compute p-value
#' @param ...  Optional arguments used to fit null model with \code{lmm.aireml} of \code{logistic.mm.aireml} function.
#' 
#' 
#' @details
#' In \code{score.variance.linear}, we consider the linear mixed model
#' \deqn{ Y = X\alpha + \gamma + \omega_1 + \ldots + \omega_k + \varepsilon }{ Y = X alpha + gamma + omega_1 + ... + omega_k + varepsilon }
#' or, in \code{score.variance.logistic}, we consider the following logistic model
#' \deqn{ \mbox{logit}(P[Y=1|X,x,\omega_1,\ldots,\omega_k]) = X\alpha + \gamma + \omega_1 + \ldots + \omega_k}{logit(P[Y=1|X,x,omega_1,...,omega_k]) = X alpha + gamma + omega_1 + ... + omega_k }
#' with \eqn{ \gamma\sim N(0,\kappa K_0)\gamma}{gamma~N(0, kappa K_0)}, \eqn{ \omega_j \sim N(0,\tau_j K_j) }{omega_j ~ N(0, tau_j K_j)},
#' \eqn{ \varepsilon \sim N(0,\sigma^2 I_n) }{epsilon ~ N(0, sigma^2 I_n)}.
#' \eqn{K_0} and \eqn{K_j} are Genetic Relationship Matrix (GRM).
#' 
#' \code{score.variance.linear } and \code{score.variance.logistic} functions permit to test
#' \deqn{H_0 : \kappa=0 \mbox{ vs } H_1 : \kappa>0}{H_0 : kappa=0 vs H_1 : kappa>0}
#' with, for linear mixed model, the score
#' \deqn{ Q = Y'P_OK_0P_0Y/2 }{ Q = Y'P_OK_0P_0Y/2 }
#' or, for logistic mixed model, the score
#' \deqn{ Q = (Y-\pi_0)'K_0(Y-\pi_0)/2 }{ Q = (Y-pi_0)'K_0(Y-pi_0)/2 }
#' where \eqn{P_0} is the last matrix \eqn{P} computed in the optimization process for null model and \eqn{\pi_0}{pi_0} the vector of fitted values under null logistic model.
#' 
#' The associated p-value is computed with Davies method.
#' 
#' 
#' In this aim, all parameters under null model are estimated with \code{lmm.aireml} or \code{logistic.mm.aireml}.
#' The p-value corresponding to the estimated score is computed using Davies method implemented in 'CompQuadForm' R package.
#' 
#' @return
#' A named list of values:
#' \item{score}{ Estimated score }
#' \item{p}{ The corresponding p-value }
#' @seealso  \code{\link{lmm.aireml}}, \code{\link{logistic.mm.aireml}}
#' 
#' @references Davies R.B. (1980) \emph{Algorithm AS 155: The Distribution of a Linear Combination of chi-2 Random Variables},
#' Journal of the Royal Statistical Society. Series C (Applied Statistics), \bold{323-333}
#' 
#' @keywords  Variance Component Test
#' @examples
#' 
#' # Load data
#' data(AGT)
#' x <- as.bed.matrix(AGT.gen, AGT.fam, AGT.bim)
#' standardize(x) <- "p"
#' 
#' # Calculate GRM et its eigen decomposition
#' K0 <- GRM(x)
#' eig <- eigen(K0)
#' eig$values <- round(eig$values, 5)
#' 
#' # generate an other positive matrix (to play the role of the second GRM)
#' set.seed(1)
#' R <- random.pm(nrow(x))
#' 
#' 
#' # simulate quantitative phenotype with two polygenic components
#' y <- lmm.simu(0.1,1,eigenK=eig)$y + lmm.simu(0.2,0,eigenK=R$eigen)$y
#' 
#' t <- score.variance.linear(K0, y, K=R$K, verbose=FALSE)
#' str(t)
#' 
#' 
#' # simulate binary phenotype with two polygenic components
#' mu <- lmm.simu(0.1,0.5,eigenK=eig)$y + lmm.simu(0.2,0,eigenK=R$eigen)$y
#' pi <- 1/(1+exp(-mu))
#' y <- 1*(runif(length(pi))<pi)
#' 
#' tt <- score.variance.logistic(K0, y, K=R$K, verbose=FALSE)
#' str(tt)
#' 
#' @export
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


