#' Duplicated SNPs
#' 
#' @description  Determines which SNPs are duplicates of previous SNPs
#' and returns their indices.
#' 
#' @param x  A bed.matrix or a data.frame
#' @param by  The criterium used to determined if SNP is duplicated.
#' 
#' 
#' @details
#' When \code{x} is a bed.matrix, the data.frame \code{x@bed} will be used.
#' The columns that will be taken in consideration
#' Are \code{id}, \code{chr}, \code{pos}, \code{A1}, and \code{A2}. Not all columns
#' are mandatory, depending on the value of \code{by}.
#' 
#' The possible values for \code{by} are \code{"chr:pos"}, \code{"chr:pos:alleles"}, \code{"id"},
#' \code{"id:chr:pos"} and \code{"id:chr:pos:alleles"}.
#' The default is \code{by = "chr:pos"}, which means that two SNPs are considered as duplicated if they have
#' same \code{chr} and \code{pos} values.
#' 
#' Currently, when using a criterium involving alleles, this function does not consider the possibility
#' of alleles swaps or reference strand flips.
#' 
#' 
#' 
#' @return An integer vector of indices of SNPs which are duplicates of previously seen SNPs.
#' @seealso  \code{\link{SNP.match}}
#' 
#' 
#' @export SNP.duplicated
SNP.duplicated <- function(x, by = "chr:pos") {
  if (is(x, "bed.matrix")) {
    x <- x@snps
  }
  if (!inherits(x, "data.frame")) {
    stop("x be a bed matrix or a data frame")
  }

  if(by == "id") {
    if(is.null(x$id)) stop("Missing id in x")
    return(.Call(`_gaston_which_duplicated_id`, PACKAGE = "gaston", x$id))
  }
  if(by == "id:chr:pos") {
    if(is.null(x$id)) stop("Missing id in x")
    if(is.null(x$chr)) stop("Missing chr in x")
    if(is.null(x$pos)) stop("Missing pos in x")
    return(.Call(`_gaston_which_duplicated_id_chr_pos`, PACKAGE = "gaston", x$id, x$chr, x$pos))
  }

  if(by == "id:chr:pos:alleles") {
    if(is.null(x$id)) stop("Missing id in x")
    if(is.null(x$chr)) stop("Missing chr in x")
    if(is.null(x$pos)) stop("Missing pos in x")
    if(is.null(x$A1) || is.null(x$A2)) stop("Missing alleles in x")
    return(.Call(`_gaston_which_duplicated_id_chr_pos_alleles`, PACKAGE = "gaston", x$id, x$chr, x$pos, x$A1, x$A2))
  }

  if(by == "chr:pos") {
    if(is.null(x$chr)) stop("Missing chr in x")
    if(is.null(x$pos)) stop("Missing pos in x")
    return(.Call(`_gaston_which_duplicated_chr_pos`, PACKAGE = "gaston", x$chr, x$pos))
  }

  if(by == "chr:pos:alleles") {
    if(is.null(x$chr)) stop("Missing chr in x")
    if(is.null(x$pos)) stop("Missing pos in x")
    if(is.null(x$A1) || is.null(x$A2)) stop("Missing alleles in x")
    return(.Call(`_gaston_which_duplicated_chr_pos_alleles`, PACKAGE = "gaston", x$chr, x$pos, x$A1, x$A2))
  }
  stop("'by' should be one of 'id', 'id:chr:pos', 'id:chr:pos:alleles', 'chr:pos', 'chr:pos:alleles'")
}
