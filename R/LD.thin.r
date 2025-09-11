#' LD thinning
#' 
#' @description  Select SNPs in LD below a given threshold.
#' 
#' @param x  A \code{\link{bed.matrix}}
#' @param threshold  The maximum LD (measured by \eqn{r^2}) between SNPs
#' @param max.dist  The maximum distance for which the LD is computed
#' @param beg  The index of the first SNP to consider
#' @param end  The index of the last SNP to consider
#' @param which.snps  Logical vector, giving which SNPs are considerd. The default is to use all SNPs
#' @param dist.unit  Distance unit in \code{max.dist}
#' @param extract  A \code{logical} indicating whether the function return a \code{bed.matrix} (\code{TRUE})
#' or a logical vector indicating which SNPs are selected (\code{FALSE})
#' @param keep  Which SNP is selected in a pair with LD above \code{threshold}
#' 
#' @details
#' The SNPs to keep are selected by a greedy algorithm. The LD is computed only for SNP pairs for which distance is inferior to
#' \code{max.dist}, expressed in number of bases if \code{dist.unit = "bases"}, in number of SNPs if \code{dist.unit = "indices"},
#' or in centiMorgan if \code{dist.unit = "cM"}.
#' 
#' The argument \code{which.snps} allows to consider only a subset of SNPs.
#' 
#' The algorithm tries to keep the largest possible number of SNPs: it is not appropriate to select tag-SNPs.
#' The algorithm takes the first SNP (snp1) on the chromosome, then goes right SNP by SNP and compute the LD.
#' When it finds a SNP (snp2) above the threshold, it haves to keep only one of (snp1, snp2).
#' The choice is made according to the parameter 'keep'.
#' \itemize{
#' \item If 'keep' is 'left' it keeps snp1 and remove snp2, and then continues to go right until max.dist is reached.
#' Once it has reached max.dist, it considers the next SNP (after snp1) that has not been removed, and iterates.
#' \item If 'keep' is 'right', it keeps snp2 and removes snp1. Then it considers the next SNP (after snp1) that has not
#' been removed, and iterates.
#' \item If 'keep' is 'random', then one of these two strategies is used, with a prob of 0.5 / 0.5.
#' }
#'
#' @return
#' If \code{extract = TRUE}, a \code{\link{bed.matrix}} extracted from \code{x} with SNPs in pairwise LD below the given threshold.
#' If \code{extract = FALSE}, a logical vector of length \code{end - beg + 1}, where \code{TRUE} indicates that
#' the corresponding SNPs is selected.
#' @seealso  \code{\link{LD}}, \code{\link{set.dist}}
#' 
#' @keywords  Linkage Disequilibrium
#' @examples
#' 
#' # Load data
#' data(TTN)
#' x <- as.bed.matrix(TTN.gen, TTN.fam, TTN.bim)
#' 
#' # Select SNPs in LD r^2 < 0.4, max.dist = 500 kb
#' y <- LD.thin(x, threshold = 0.4, max.dist = 500e3)
#' y
#' 
#' # Verifies that there is no SNP pair with LD r^2 > 0.4
#' # (note that the matrix ld.y has ones on the diagonal)
#' ld.y <- LD( y, lim = c(1, ncol(y)) )
#' sum( ld.y > 0.4 )  
#' 
#' @export LD.thin
LD.thin <- function(x, threshold, max.dist = 500e3, beg = 1L, end = ncol(x), which.snps,
                    dist.unit = c("bases", "indices", "cM"), extract = TRUE, keep = c("left", "right", "random")) {

  if(missing(which.snps)) which.snps <- rep(TRUE, end-beg+1)

  if(!is.logical(which.snps) | length(which.snps) != end-beg+1)
    stop("which.snps must be a Logical vector of length end - beg + 1")

  if(is.null(x@mu) | is.null(x@sigma))
    stop("LD.thin needs mu and sigma to be set for LD computation (use set.stats)")

  # ne pas considérer les SNPs monomorphes ou qui ont un callrate nul
  I <- beg:end
  which.snps <- which.snps & (x@snps$callrate[I] > 0) & (x@snps$maf[I] > 0)

  dist.unit <- match.arg(dist.unit)
  if(dist.unit == "bases") {
    pos <- as.integer(x@snps$pos)
  } else if(dist.unit == "indices") {
    pos <- seq_len(ncol(x))
  } else { # cM = conversion en micro morgans !
    check.dist(x) # check if 'dist' is set
    if(max.dist > 100)
      warning("max.dist value seems very high for dist.unit = \"cM\"")
    pos <- as.integer( x@snps$dist * 1e4 )
    max.dist <- as.integer(max.dist*1e4)
  }
  if( all(x@snps$pos == x@snps$pos[1]) )
    stop("Position of SNPs must be available")

  keep <- match.arg(keep)
  if(keep == "left") {
    w <- .Call(`_gaston_ld_thin_left`, x@bed, x@mu, x@sigma, threshold, pos, 
          as.integer(x@snps$chr), as.integer(max.dist), as.integer(beg)-1L, as.integer(end)-1L,
          which.snps)
  } else if (keep == "right"){
    w <- .Call(`_gaston_ld_thin_right`, x@bed, x@mu, x@sigma, threshold, pos, 
          as.integer(x@snps$chr), as.integer(max.dist), as.integer(beg)-1L, as.integer(end)-1L,
          which.snps)
  } else {
    w <- .Call(`_gaston_ld_thin_random`, x@bed, x@mu, x@sigma, threshold, pos, 
          as.integer(x@snps$chr), as.integer(max.dist), as.integer(beg)-1L, as.integer(end)-1L,
          which.snps)
  }

  if(!extract) return(w)
  x[ , seq(beg,end)[w] ]
}

