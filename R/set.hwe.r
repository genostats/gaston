#' Hardy-Weinberg Equilibrium
#' 
#' @description
#' Returns an updated \code{\link{bed.matrix}} with a new variable for the \eqn{p}-values of an
#' Hardy-Weinberg Equilibrium test.
#' 
#' @param x  A \code{\link{bed.matrix}}
#' @param method  The method to use, either "chisquare" or "exact"
#' @param verbose  If \code{TRUE}, display information on the function actions
#' 
#' @details
#' Two tests of Hardy-Weinberg Equilibrium are proposed:
#' \itemize{
#' \item if \code{method = "chisquare"}, the good old Chi-square test
#' \item if \code{method = "exact"}, Haldane's exact test (see Wigginton et al)
#' }
#' 
#' The function \code{set.stats} will be called first if necessary.
#' 
#' The \eqn{p}-value is set to \eqn{1.0} for SNPs on chromosomes Y and MT. For SNPs on
#' chromosomes X, currently, the test is performed using only the genotypic counts of women.
#' 
#' @return
#' A \code{\link{bed.matrix}} similar to \code{x}, with a new variable \code{x@snps$hwe}
#' containing the \eqn{p}-values for each SNP.
#' @seealso  \code{\link{set.stats}}, \code{\link{set.genomic.sex}}
#' 
#' @references  Wigginton, J. E., Cutler, D. J., & Abecasis, G. R. (2005). \emph{A note on exact tests of Hardy-Weinberg equilibrium}. The American Journal of Human Genetics, \bold{76(5)}, \bold{887-893}
#' 
#' @keywords  Hardy-Weinberg   P-value
#' @examples
#' 
#' # Load data
#' data(LCT)
#' x <- as.bed.matrix(LCT.gen, LCT.fam, LCT.bim)
#' 
#' # Compute Hardy-Weinberg p-values
#' x <- set.hwe(x)
#' head( x@snps[,c("id","hwe")] )
#' 
#' @export set.hwe
set.hwe <- function(x, method = c("chisquare", "exact"), verbose = getOption("gaston.verbose",TRUE)) {
  if( !all(c("N0", "N1", "N2") %in% names(x@snps) )) {
    if(verbose) cat("Computing basic stats\n")
    x <- set.stats(x)
  }

  w.a <- is.autosome(x@snps$chr)
  w.x <- is.chr.x(x@snps$chr)

  hwe <- rep(1.0, ncol(x) )
  method <- match.arg(method)
  if(method == 'chisquare') {
    if(verbose) cat("Computing HW chi-square p-values\n")
    hwe[w.a] <- .Call(`_gaston_hwe_chi`, PACKAGE = "gaston", x@snps$N0[w.a], x@snps$N1[w.a], x@snps$N2[w.a])
    hwe[w.x] <- .Call(`_gaston_hwe_chi`, PACKAGE = "gaston", x@snps$N0.f[w.x], x@snps$N1.f[w.x], x@snps$N2.f[w.x])
  } else {
    if(verbose) cat("Computing HW exact test p-values\n")
    hwe[w.a] <- .Call(`_gaston_hwe`, PACKAGE = "gaston", x@snps$N0[w.a], x@snps$N1[w.a], x@snps$N2[w.a])
    hwe[w.x] <- .Call(`_gaston_hwe`, PACKAGE = "gaston", x@snps$N0.f[w.x], x@snps$N1.f[w.x], x@snps$N2.f[w.x])
  }
  x@snps$hwe <- hwe
  x
}



