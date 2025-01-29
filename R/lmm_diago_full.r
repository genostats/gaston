lmm.diago.profile.likelihood <- function(tau, s2, h2, Y, X = matrix(1, nrow=length(Y)), eigenK, p = 0) {
  if( any(is.na(Y)) ) 
    stop('Missing data in Y.')

  if(!is.null(X)) {
    if( length(Y)<(ncol(X)+max(p)) ) 
      stop('The total number of covariables and PCs as fixed effects cannot exceed the number of observations.')
    if( length(Y)!=nrow(X) ) 
      stop('Length of Y and the number of rows of X differ.')
  } else if( length(Y) < max(p) )
    stop('The total number of covariables and PCs as fixed effects cannot exceed the number of observations.')

  Sigma <- eigenK$values
  if( length(Y)!=length(Sigma) )
    stop('Length of Y and number of eigenvalues differ.')

  if(!missing(tau) & !missing(s2)) 
    if(is.null(X))
      .Call(`_gaston_diago_full_likelihood2_nocovar`, PACKAGE = "gaston", tau, s2, p, Y, Sigma, eigenK$vectors)
    else
      .Call(`_gaston_diago_full_likelihood2`, PACKAGE = "gaston", tau, s2, p, Y, X, Sigma, eigenK$vectors)
  else if(!missing(h2))
    if(is.null(X))
      .Call(`_gaston_diago_full_likelihood1_nocovar`, PACKAGE = "gaston", h2, p, Y, Sigma, eigenK$vectors)
    else
      .Call(`_gaston_diago_full_likelihood1`, PACKAGE = "gaston", h2, p, Y, X, Sigma, eigenK$vectors)
  else
    stop("To compute likelihood, provide tau and s2 values, or h2 value")
}


