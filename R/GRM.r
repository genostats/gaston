#' Genetic Relationship Matrix
#' 
#' @description
#' Compute the Genetic Relationship Matrix
#' 
#' @param x  A \code{\link{bed.matrix}}
#' @param which.snps  Logical vector, giving which snps to use in the computation.  The default is to use all autosomal SNPs
#' @param autosome.only  If \code{TRUE}, only autosomal SNPs will be considered.
#' @param chunk  Parameter for the parallelization: how many SNPs are treated by each task
#' 
#' @details
#' The Genetic Relationship Matrix (GRM) is computed by the formula \eqn{{1\over q}XX'}{XX'/q},
#' with \eqn{X} the standardized genotype matrix and \eqn{q} the number of SNPs
#' (\code{ncol(x)}).
#' 
#' If \code{x} is not standardized before this computation, the function
#' will use \code{standardize(x) <- "p"} by default.
#' 
#' @return  The GRM is a symmetric square matrix of dimension equal to the number of individuals.
#' Each entry can be interpreted as an estimated kinship coefficient between individuals, although some
#' authors might disagree. Note in particular that some entries will be negative.
#' @seealso  \code{\link{DM}}, \code{\link{reshape.GRM}}, \code{\link{lmm.aireml}}, \code{\link{lmm.diago}}, \code{\link{standardize}}, \code{\link{bed.loadings}}
#' 
#' @keywords  Genetic Relationship Matrix
#' @examples
#' 
#' # load chr2 data set (~10k SNPs in low LD)
#' x <- read.bed.matrix( system.file("extdata", "chr2.bed", package="gaston") )
#' 
#' # Compute Genetic Relationship Matrix
#' K <- GRM(x)
#' dim(K)
#' 
#' @export GRM
GRM <- function(x, which.snps, autosome.only = TRUE, chunk = 1L) {
  if(missing(which.snps)) which.snps <- rep(TRUE, ncol(x))
  if(autosome.only) 
    which.snps <- which.snps & is.autosome(x@snps$chr)

  if(!is.logical(which.snps) | length(which.snps) != ncol(x))
    stop("which.snps must be a Logical vector of length ncol(x)")
  
  if(!x@standardize_mu_sigma & !x@standardize_p) {
    if(!is.null(x@p)) x@standardize_p <- TRUE
    else stop("Can't standardize x for GRM computation (use set.stat)\n")
  }

  if(x@standardize_mu_sigma) {
    w <- ifelse(x@sigma == 0, 0, 1/x@sigma/sqrt(sum(which.snps)-1))   ### BEWARE q-1 !!!
    K <- .Call(`_gaston_Kinship_w`, PACKAGE = "gaston", x@bed, x@mu[which.snps], w[which.snps], which.snps, chunk) 
  } else { 
    K <- .Call(`_gaston_Kinship_pw`, PACKAGE = "gaston", x@bed, x@p[which.snps], which.snps, FALSE, chunk)
  }

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
