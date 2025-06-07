#' LD clumping
#' 
#' @description  Construct group of SNPs in LD with 'top associated SNPs'
#' 
#' @param x  A \code{\link{bed.matrix}}
#' @param p  A vector of p-values, or a data frame including p-values, such as sent back by \code{\link{association.test}}
#' @param r2.threshold  The maximum LD (measured by \eqn{r^2}) between SNPs in a group
#' @param p.threshold  The threshold used to define associated SNPs
#' @param max.dist  The maximum distance for which the LD is computed
#' 
#' @details
#' The p-values provided through argument \code{p} are assumed to correspond to the result of an association test with the SNPs of \code{x}.
#' 
#' The aim of the function is to construct cluster of SNPs in strong LD with associated SNPs.
#' The algorithm first seeks the SNP with the lowest p-value (below \code{p.threshold}) ; this SNP will be the 'index' of a cluster.
#' The corresponding cluster is constructed by aggregating SNPs that are in LD (above \code{r2.threshold}) with the index. The cluster's name
#' is the position of the index SNP.
#' The processus is repeated on the SNPs which are not yet attributed to a cluster, until there is no associated SNP
#' (ie SNP with a p-value below \code{threshold}) left.
#' The remaining SNPs are attributed to cluster 0.
#' 
#' The LD is computed only for SNP pairs for which distance is inferior to \code{max.dist}, expressed in number of bases: above this
#' distance it is assumed to be null.
#' 
#' @return
#' If \code{p} was a data frame, then the function returns the same data frame with to extra columns, \code{cluster} and \code{is.index}.
#' If \code{p} was a vector of p-values, it returns a data frame with columns \code{chr}, \code{id}, \code{pos}, \code{p}, \code{cluster}
#' and \code{is.index}.
#' @seealso  \code{\link{LD}}, \code{\link{LD.thin}}
#' 
#' @keywords  Linkage Disequilibrium
#' @examples
#' 
#' # Construct a bed matrix
#' x <- as.bed.matrix(TTN.gen, TTN.fam, TTN.bim)
#' standardize(x) <- "p"
#'      
#' # simulate quantitative phenotype with effect of SNPs #108 and #631
#' beta <- numeric(ncol(x))
#' beta[c(108,631)] <- 0.5
#' set.seed(1)
#' y <- x %*% beta + rnorm(nrow(x))
#'      
#' # association test with linear model 
#' test <- association.test(x, y, method="lm", response = "quanti")
#' 
#' # LD clumping
#' test <- LD.clump(x, test, r2.threshold = 0.25, p.threshold = 1e-8)
#' 
#' # use as.factor for a quick-and-dirty cluster colouring on the manhattan plot 
#' manhattan(test, col = as.factor(test$cluster), pch = 20)
#' 
#' @export LD.clump
LD.clump <- function(x, p, r2.threshold, p.threshold, max.dist = 500e3) {
  if(!is(x, "bed.matrix")) 
    stop("x is not a bed matrix")
  if(is.data.frame(p)) {
    if(!all(p$chr == x@ped$chr, na.rm = TRUE) | !all(p$pos == x@ped$pos, na.rm = TRUE))
      stop("Unmatching SNPs")
    a <- p
    p <- p$p
  } else {
    a <- NULL
  }
  if(length(p) != ncol(x))
    stop("Dimensions mismatch")

  or <- order(p) - 1L # -1 pour des indices qui démarrent à 0 ds la fction c++
  if(!missing(p.threshold)) {
    or <- or[ p[or + 1L] < p.threshold ]
    if(length(or) == 0) stop("too stringent p.threshold ?")
  }

  I <- .Call(`_gaston_ld_clump`, PACKAGE = "gaston", x@bed, x@mu, x@sigma, r2.threshold, x@snps$pos, x@snps$chr, max.dist, or);

  if(is.null(a))
    a <- data.frame( chr = x@snps$chr, id = x@snps$id, pos = x@snps$pos, p = p, cluster = I+1)
  else {
    a$cluster <- I+1
  }
  a$is.index <- (a$cluster == 1:nrow(a))
  a
}


