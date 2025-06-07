#' Dominance Matrix
#' 
#' @description
#' Compute the Dominance Matrix
#' 
#' @param x  A \code{\link{bed.matrix}}
#' @param which.snps  Logical vector, giving which snps to use in the computation.  The default is to use all autosomal SNPs
#' @param autosome.only  If \code{TRUE}, only autosomal SNPs will be considered.
#' @param chunk  Parameter for the parallelization: how many SNPs are treated by each task
#' 
#' 
#' @details
#' The Dominance Matrix (DM) gives for each pair of individuals an estimation of
#' their probability of sharing two alleles Identical By Descent.
#' 
#' It is computed by a moment estimator,
#' \eqn{{1\over q} ZZ'}{ZZ'/q} with \eqn{Z} the matrix with entries
#' \eqn{p \over 1-p}{p/(1-p)}, \eqn{-1}, \eqn{1-p \over p}{(1-p)/p} according to the
#' values 0, 1, 2 in the genotype matrix \code{x} (here \eqn{p} is the
#' frequency of the alternate allele, and \eqn{q} is the number of SNPs
#' (\code{ncol(x)}).
#' 
#' @return  A symmetric square matrix of dimension equal to the number of individuals.
#' Each entry can be interpreted as an estimated probability of sharing two alleles IBD
#' — as it is a moment estimator, the value can (and will) fall outside of the range (0,1).
#' @seealso  \code{\link{GRM}}, \code{\link{reshape.GRM}}
#' 
#' @keywords  Dominance Matrix
#' @examples
#' 
#' # load chr2 data set (~10k SNPs in low LD)
#' x <- read.bed.matrix( system.file("extdata", "chr2.bed", package="gaston") )
#' 
#' # Compute Dominance Matrix
#' D <- DM(x)
#' dim(D)
#' 
#' 
#' @export DM
DM <- function(x, which.snps, autosome.only = TRUE, chunk = 1L) {
  if(missing(which.snps)) which.snps <- rep(TRUE, ncol(x))
  if(autosome.only) 
    which.snps <- which.snps & is.autosome(x@snps$chr)

  if(!is.logical(which.snps) | length(which.snps) != ncol(x))
    stop("which.snps must be a Logical vector of length ncol(x)")
  
  if(x@standardize_mu_sigma)
    warning("For Dominance Matrix, p standardization is used\n")
  if(is.null(x@p))
    stop("Can't standardize x for DM computation (use set.stat)\n")

  K <- .Call(`_gaston_Kinship_pw`, PACKAGE = "gaston", x@bed, x@p[which.snps], which.snps, TRUE, chunk)

  if(!is.null(x@ped$id)) {
    if(anyDuplicated(x@ped$id) == 0) 
      rownames(K) <- colnames(K) <- x@ped$id
    else {
      nn <- paste(x@ped$famid, x@ped$id, sep = ":")
      if(anyDuplicated(nn) == 0)
        rownames(K) <- colnames(K) <- nn
    }
  }

  K
}

