# slot p
setGeneric('p', function(x) standardGeneric('p'), package="gaston")

#' @exportMethod p
setMethod('p', signature =  'bed.matrix', def = function(x) x@p)

#' @exportMethod "p<-"
setGeneric('p<-', function(x,value) standardGeneric('p<-'), package="gaston")
setReplaceMethod('p', 'bed.matrix', function(x,value) {
  if(length(value) != ncol(x)) stop("dimensions mismatch")
  x@p <- value; 
  x} )

# slot mu
setGeneric('mu', function(x) standardGeneric('mu'), package="gaston")

#' @exportMethod mu
setMethod('mu', signature =  'bed.matrix', def = function(x) x@mu)

setGeneric('mu<-', function(x,value) standardGeneric('mu<-'), package="gaston")

#' @exportMethod "mu<-"
setReplaceMethod('mu', 'bed.matrix', function(x,value) {
  if(length(value) != ncol(x)) stop("dimensions mismatch")
  x@mu <- value; 
  x} )

# slot sigma
setGeneric('sigma', function(x) standardGeneric('sigma'), package="gaston")

#' @exportMethod sigma
setMethod('sigma', signature =  'bed.matrix', def = function(x) x@sigma)

#' @exportMethod "sigma<-"
setGeneric('sigma<-', function(x,value) standardGeneric('sigma<-'), package="gaston")
setReplaceMethod('sigma', 'bed.matrix', function(x,value) {
  if(length(sigma) != ncol(x)) stop("dimensions mismatch")
  x@sigma <- value; 
  x} )

# center_scale
setGeneric('standardize', function(x) standardGeneric('standardize'), package="gaston")

#' @exportMethod standardize
setMethod('standardize', signature =  'bed.matrix', 
  def = function(x) if(x@standardize_p) "p" else if(x@standardize_mu_sigma) "mu_sigma" else "none"
)

#' @exportMethod "standardize<-"
setGeneric('standardize<-', function(x,value) standardGeneric('standardize<-'), package="gaston")
setReplaceMethod('standardize', 'bed.matrix', function(x,value=c("p","mu_sigma","none")) { 
  switch(match.arg(value), 
    p        = { if(length(x@p) != ncol(x)) 
                   stop("x@p must be defined and have the appropriate length");  
                 x@standardize_p <- TRUE;  x@standardize_mu_sigma <- FALSE;
               } ,
    mu_sigma = { if(length(x@mu) != ncol(x) | length(x@sigma) != ncol(x)) 
                   stop("x@mu and x@sigma must be defined and have the appropriate length"); 
                 x@standardize_p <- FALSE; x@standardize_mu_sigma <- TRUE; 
               } ,
    none     = {x@standardize_p <- FALSE; x@standardize_mu_sigma <- FALSE;} 
  );
  x } )


