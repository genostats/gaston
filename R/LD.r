#' Linkage Disequilibrium
#' 
#' @description  Compute Linkage Disequilibrium (LD) between given SNPs.
#' 
#' @param x  A \code{\link{bed.matrix}}
#' @param lim  Range of SNPs for which the LD is computed
#' @param lim2  (Optional) Second range of SNPs (see Details)
#' @param measure  The LD measure
#' @param trim  \code{Logical}. If \code{TRUE}, the values above 1 or below -1 are replaced by 1 and -1 respectively.
#' 
#' @details
#' If \code{lim2} is missing, the LD is computed between all SNPs with indices between \code{lim[1]} and \code{lim[2]};
#' else, the LD is computed between the SNPs in the range given by \code{lim} and those in the range given by \code{lim2}.
#' 
#' Note that the LD estimates are moment estimates (which are less precise than Maximum Likelihood Estimates).
#' If \code{standardize(x) = "none"}, \code{x} will be standardized
#' using \code{x@mu} and \code{x@sigma}. If \code{standardize(x) = "p"}, the moment estimates can produce \eqn{r}
#' values outside of the range \eqn{[-1;1]}, hence the parameter \code{trim}. We recommend to set
#' \code{standardize(x) <- "mu"} (trimming can still be necessary due to rounding errors).
#' 
#' @return
#' A matrix of LD values.
#' @seealso  \code{\link{LD.thin}},  \code{\link{LD.plot}}
#' 
#' @keywords  Linkage Disequilibrium
#' @examples
#' 
#' # Load data
#' data(AGT)
#' x <- as.bed.matrix(AGT.gen, AGT.fam, AGT.bim)
#' 
#' # Compute LD
#' ld.x <- LD(x, c(1,ncol(x)))
#' 
#' # Plot a tiny part of the LD matrix
#' LD.plot( ld.x[1:20,1:20], snp.positions = x@snps$pos[1:20] )
#' 
#' @export LD
LD <- function(x, lim, lim2, measure = c("r2", "r", "D"), trim = TRUE) {

  if(!is(x,"bed.matrix")) stop('x must be a bed matrix')
  if(!is.vector(lim))  stop('lim must be a vector of length 2')
  if(length(lim) == 1) lim = c(lim, lim)
  if(length(lim) != 2) stop('lim must be a vector of length 2')
  if(any(lim < 1) | any(lim > ncol(x)) | lim[1] > lim[2]) stop("inappropriate range in lim")

  if(!missing(lim2)) {
    if(!is.vector(lim2))  stop('lim2 must be a vector of length 2')
    if(length(lim2) == 1) lim2 = c(lim2, lim2)
    if(length(lim2) != 2) stop('lim2 must be a vector of length 2')
    if(any(lim2 < 1) | any(lim2 > ncol(x)) | lim2[1] > lim2[2]) stop("inappropriate range in lim2")
  }

  measure <- match.arg(measure) 

  if(!x@standardize_mu_sigma & !x@standardize_p) {
     if(!is.null(x@mu) & !is.null(x@sigma))
       x@standardize_mu_sigma <- TRUE
     else
       stop("Can't standardize x for LD computation (use set.stats(x))\n")
  }

  if(measure == "D") {
    x@standardize_mu_sigma <- TRUE
    x@sigma <- rep(sqrt(2), ncol(x))
  }  

  if(missing(lim2) || all(lim2 == lim)) {  
    if(x@standardize_mu_sigma) 
      a <- .Call(`_gaston_LD_`, PACKAGE = "gaston", x@bed, x@mu, x@sigma, lim[1]-1, lim[2]-1)
    else if(x@standardize_p) {
      warning("Moment estimates of LD using a p-standardized matrix can be outside of the range [-1,1]")
      a <- .Call(`_gaston_LD_p`, PACKAGE = "gaston", x@bed, x@p, lim[1]-1, lim[2]-1)
    }
    rownames(a) <- colnames(a) <- x@snps$id[seq(lim[1], lim[2])]
  } else { 
    if(x@standardize_mu_sigma) 
      a <- .Call(`_gaston_LD_chunk`, PACKAGE = "gaston", x@bed, x@mu, x@sigma, lim[1]-1, lim[2]-1, lim2[1]-1, lim2[2]-1)
    else if(x@standardize_p)
      a <- .Call(`_gaston_LD_chunk_p`, PACKAGE = "gaston", x@bed, x@p, lim[1]-1, lim[2]-1, lim2[1]-1, lim2[2]-1)
    rownames(a) <- x@snps$id[seq(lim[1], lim[2])]
    colnames(a) <- x@snps$id[seq(lim2[1], lim2[2])]
  }

  if(measure != "D" & trim) {
    a[a < -1] <- -1; 
    a[a >  1] <-  1;
  }

  if(measure == "r2") a <- a**2

  a
}

