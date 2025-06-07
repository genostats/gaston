#' SNP loadings
#' 
#' @description  Compute the loadings corresponding to given PCs.
#' 
#' @param x  A \code{\link{bed.matrix}}
#' @param pc  A matrix with Principal Components in columns
#' 
#' @return  A matrix with the corresponding loadings in columns.
#' 
#' 
#' @keywords  loadings
#' @examples
#' 
#' # load chr2 data set (~10k SNPs in low LD)
#' x <- read.bed.matrix( system.file("extdata", "chr2.bed", package="gaston") )
#' 
#' # Compute Genetic Relationship Matrix
#' standardize(x) <- "p"
#' K <- GRM(x)
#' 
#' # Eigen decomposition
#' eiK <- eigen(K)
#' # deal with small negative eigen values
#' eiK$values[ eiK$values < 0 ] <- 0
#' 
#' # Note: the eigenvectors are normalized, to compute 'true' PCs
#' # multiply them by the square root of the associated eigenvalues
#' PC <- sweep(eiK$vectors, 2, sqrt(eiK$values), "*")
#' 
#' # Compute loadings for the 2 first PCs 
#' # one can use PC[,1:2] instead of eiK$vectors[,1:2] as well
#' L <- bed.loadings(x, eiK$vectors[,1:2])
#' dim(L)
#' head(L)
#' 
#' # the loadings are normalized
#' colSums(L**2)
#' 
#' # Verify that these are loadings
#' head( (x %*% L) / sqrt(ncol(x)-1) )
#' head( PC[,1:2] )
#' 
#' 
#' @export bed.loadings
bed.loadings <- function(x, pc) {
  if(is.vector(pc)) dim(pc) <- c(length(pc),1)
  if(!is.matrix(pc)) stop("pc must be a matrix of PCs")
  if(!is(x, "bed.matrix")) stop("x must be a bed.matrix")
  if(!x@standardize_mu_sigma & !x@standardize_p) {
    if(!is.null(x@p)) 
      x@standardize_p <- TRUE
    else stop("Can't standardize x for LD computation (use set.stat)\n")
 }

  if(x@standardize_mu_sigma)
    a <- .Call(`_gaston_m4_pc_to_loading_ms`, PACKAGE = "gaston", x@bed, x@mu, x@sigma, pc)
  else if(x@standardize_p)
    a <- .Call(`_gaston_m4_pc_to_loading_p`, PACKAGE = "gaston", x@bed, x@p, pc)

  a <- a/sqrt(ncol(x))
  a <- sweep(a, 2L, sqrt(colSums(a**2)), "/")
  if(!is.null(x@snps$id)) {
    if(anyDuplicated(x@snps$id) == 0)
      rownames(a) <- x@snps$id
  }
  a
}

