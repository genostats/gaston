#' Evaluation of a condition on SNPS or individuals in a \code{\link{bed.matrix}}
#' 
#' @description
#' Evaluate a condition and return logical vector or indices
#' 
#' @param x  A \code{\link{bed.matrix}}
#' @param condition  Condition used to select SNPs
#' @param na.to.false  If \code{TRUE}, \code{NA}s are replaced by \code{FALSE}
#' 
#' 
#' @details
#' The conditions can involve global variables and all variables defined
#' in the data frame \code{x@snps}, in particular for \code{test.snps} and \code{which.snps}
#' \itemize{
#' \item \code{chr}, \code{id}, \code{dist}, \code{pos}, \code{A1}, \code{A2}
#' \item If basic stats have been computed (see \code{\link{set.stats}}), \code{N0}, \code{N1}, \code{N2}, \code{NAs}, \code{callrate}, \code{maf}, \code{hz}, etc.
#' \item If Hardy-Weinberg Equilibrium test has been performed (see \code{\link{set.hwe}}), \code{hwe}.
#' }
#' and for  \code{test.inds} and \code{which.inds}
#' \itemize{
#' \item \code{famid}, \code{id}, \code{father}, \code{mother}, \code{sex}, \code{pheno}
#' \item If basic stats have been computed (see \code{\link{set.stats}}),
#' \code{N0}, \code{N1}, \code{N2}, \code{NAs}, \code{callrate}, etc.
#' }
#' 
#' 
#' 
#' @return  \code{test.snps} and \code{test.inds} return a logical vector of length \code{ncol(x)} and \code{nrow(x)} respectively. \code{which.snps(x, condition)} is
#' equivalent to \code{which(test.snps(x, condition))} and \code{which.inds(x, condition)} to \code{which(test.inds(x, condition))}.
#' @seealso  \code{\link{select.snps}}, \code{\link{select.inds}}, \code{\link{set.stats}}, \code{\link{set.hwe}}
#' 
#' 
#' @examples
#' 
#' # Load data
#' data(LCT)
#' x <- as.bed.matrix(LCT.gen, LCT.fam, LCT.bim)
#' 
#' # SNPs and individuals with a callrate < 100%
#' w <- test.snps(x, callrate < 1)
#' table(w)
#' which.snps(x, callrate < 1)
#' which.inds(x, callrate < 1)
#' 
#' @export
test.snps <- function(x, condition, na.to.false = TRUE) {
  if(!is(x, "bed.matrix")) stop("x is not a bed.matrix")
  w <- eval(substitute(condition), x@snps, parent.frame())
  miss <- is.na(w)
  if(na.to.false & sum(miss)>0)
    return(w & !miss)
  w
}

#' @rdname test.snps
#' @export
test.inds <- function(x, condition, na.to.false = TRUE) {
  if(!is(x, "bed.matrix")) stop("x is not a bed.matrix")
  w <- eval(substitute(condition), x@ped, parent.frame())
  miss <- is.na(w)
  if(sum(miss)>0)
    return(w & !miss)
  w
}

#' @rdname test.snps
#' @export
which.snps <- function(x, condition) {
  if(!is(x, "bed.matrix")) stop("x is not a bed.matrix")
  which(eval(substitute(condition), x@snps, parent.frame()))
}

#' @rdname test.snps
#' @export
which.inds <- function(x, condition) {
  if(!is(x, "bed.matrix")) stop("x is not a bed.matrix")
  which(eval(substitute(condition), x@ped, parent.frame()))
}


