#' Subsetting from a \code{\link{bed.matrix}}
#' 
#' @description
#' Returns subset of individuals satisfying a condition.
#' 
#' @param x  A \code{\link{bed.matrix}}
#' @param condition  Condition used to select individuals
#' 
#' 
#' @details
#' The conditions can involve global variables and all variables defined
#' in the data frame \code{x@ped}, in particular
#' \itemize{
#' \item \code{famid}, \code{id}, \code{father}, \code{mother}, \code{sex}, \code{pheno}
#' \item If basic stats have been computed (see \code{\link{set.stats}}),
#' \code{N0}, \code{N1}, \code{N2}, \code{NAs}, \code{callrate}, etc.
#' }
#' If some condition evaluate to \code{NA} (e.g. \code{sex == 1} when \code{sex} is undefined for some individuals),
#' a warning is issued and the corresponding individuals are removed.
#' 
#' 
#' @return  A \code{\link{bed.matrix}} similar to \code{x}, containing the selected individuals only
#' @seealso  \code{\link{select.snps}}, \code{\link{set.stats}}
#' 
#' 
#' @examples
#' 
#' # Load data
#' data(LCT)
#' x <- as.bed.matrix(LCT.gen, LCT.fam, LCT.bim)
#' 
#' # Select individuals with a call rate > 95% 
#' # and more than 5% of heterozygous genotypes
#' y <- select.inds(x, callrate > 0.95 & N1/(N0+N1+N2) > 0.05)
#' y
#' 
#' @export select.inds
select.inds <- function(x, condition) {
  if(!is(x, "bed.matrix")) stop("x is not a bed.matrix")
  w <- eval(substitute(condition), x@ped, parent.frame())
  miss <- is.na(w)
  if(sum(miss)>0) {
    warning(paste(sum(miss), 'individual(s) with undefined condition are removed from bed.matrix'))
    w <- w & !miss
  }
  x[w,]
}

