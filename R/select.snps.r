#' Subsetting from a \code{\link{bed.matrix}}
#' 
#' @description
#' Returns subset of SNPs satisfying a condition.
#' 
#' @param x  A \code{\link{bed.matrix}}
#' @param condition  Condition used to select SNPs
#' 
#' @details
#' The conditions can involve global variables and all variables defined
#' in the data frame \code{x@snps}, in particular
#' \itemize{
#' \item \code{chr}, \code{id}, \code{dist}, \code{pos}, \code{A1}, \code{A2}
#' \item If basic stats have been computed (see \code{\link{set.stats}}), 
#'       \code{N0}, \code{N1}, \code{N2}, \code{NAs}, \code{callrate}, \code{maf}, \code{hz}, etc.
#' \item If Hardy-Weinberg Equilibrium test has been performed (see \code{\link{set.hwe}}), \code{hwe}.
#' }
#' 
#' If some condition evaluate to \code{NA} (e.g. \code{maf > 0} when \code{maf} is undefined for some SNPs),
#' a warning is issued and the corresponding SNPs are removed.
#' 
#' @return  A \code{\link{bed.matrix}} similar to \code{x}, containing the selected SNPs only
#' @seealso  \code{\link{select.snps}}, \code{\link{set.stats}}, \code{\link{set.hwe}}
#' 
#' @examples
#' 
#' # Load data
#' data(LCT)
#' x <- as.bed.matrix(LCT.gen, LCT.fam, LCT.bim)
#' 
#' # Select SNPs with a maf > 5%
#' y <- select.snps(x, maf > 0.05)
#' y
#' 
#' @export select.snps
select.snps <- function(x, condition) {
  if(!is(x, "bed.matrix")) stop("x is not a bed.matrix")
  w <- eval(substitute(condition), x@snps, parent.frame())
  miss <- is.na(w)
  if(sum(miss)>0) {
    warning(paste(sum(miss), 'SNP(s) with undefined condition are removed from bed.matrix'))
    w <- w & !miss
  }
  x[,w]
}

