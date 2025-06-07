#' Remove duplicated SNPs
#' 
#' @description  Remove duplicated SNPs, taking into account possible genotype mismatches
#' 
#' @param x  A bed.matrix
#' @param by  The criterium used to determine duplicates
#' @param na.keep  If \code{TRUE}, duplicated genotypes which are missing for at
#' least one SNP are set to \code{NA}.
#' @param incomp.rm  If \code{TRUE}, duplicated SNPs with allele incompatibility are
#' removed.
#' 
#' 
#' @details
#' Positions of duplicated SNPs are determined using \code{\link{SNP.duplicated}}
#' using parameter \code{by} (we recommend to use \code{"chr:pos"}, the default).
#' 
#' Then the function considers the possibility of alleles swaps or reference strand flips.
#' In case of allele incompatibility, the SNPs can be removed or not (according to \code{incomp.rm}
#' parameter).
#' 
#' When alleles can be matched, only one of the two SNPs is conserved. If there are
#' genotype incompatibilities between the duplicates for some individuals, these genotypes are set
#' to \code{NA}. The parameter \code{na.keep} settles the case of genotypes missing in one
#' of the SNPs.
#' 
#' Moreover the function takes special care of SNP with possible alleles \code{"0"}.
#' This case occurs for monomorphic SNPs, when data are read from a \code{.ped} file; for
#' example, a whole column of \code{A A}'s will result in a SNP with alleles \code{"A"} and
#' \code{"0"}. If there's a duplicate of the SNP with a few, says, \code{A C}'s in it,
#' it will have alleles \code{"A"} and \code{"C"}. In that case, \code{\link{SNP.duplicated}}
#' with \code{by = "chr:pos:alleles"} will not consider these SNPs as duplicates.
#' 
#' @return A bed.matrix without duplicated SNPs.
#' @seealso  \code{\link{SNP.match}}, \code{\link{SNP.duplicated}}, \code{\link{dupli}}
#' 
#' 
#' @examples
#' 
#' # Use example data of 10 individuals with 7 duplicated SNPs
#' data(dupli)
#' x <- as.bed.matrix(dupli.gen, fam = dupli.ped, bim = dupli.bim)
#' 
#' # There are any duplicated positions:
#' dupli.bim
#' 
#' x1 <- SNP.rm.duplicates(x)
#' # By default (na.keep = TRUE), as soon as the genotype is missing
#' # in one of the SNPs it is set to missing 
#' # (here looking at duplicated SNPs 2a and 2b)
#' as.matrix(x[,2:3])
#' as.matrix(x1[,2])
#' 
#' # With na.keep = FALSE 
#' x2 <- SNP.rm.duplicates(x, na.keep = FALSE)
#' as.matrix(x2[,2])
#' 
#' # Let's examinate SNP 3.a and 3.b (swapped alleles)
#' as.matrix(x[,4:5])
#' as.matrix(x1[,3])
#' as.matrix(x2[,3])
#' 
#' # and so on... (see also ?dupli)
#' 
#' @export SNP.rm.duplicates
SNP.rm.duplicates <- function(x, by = "chr:pos", na.keep = TRUE, incomp.rm = TRUE) {
  if (!is(x, "bed.matrix")) stop("x should be a bed.matrix")
  
  b <- strsplit(by,':')[[1]]
  if ('alleles' %in% b) b <- c(b[b!='alleles'], 'A1', 'A2')
  
  # where are duplicated SNPs
  dupli <- SNP.match(x@snps[,b,drop=FALSE], x@snps[SNP.duplicated(x, by=by),], by = by)$index
	
  a <- .Call(`_gaston_alleles_duplicated`, PACKAGE = "gaston", x@snps, dupli)
  
  if (incomp.rm) w <- which(a$keep == TRUE)
  else  w <- which(is.na(a$keep) | a$keep == TRUE)
  
  bed <- .Call(`_gaston_duplicated_remove`, PACKAGE = "gaston", x@bed, dupli, a$keep, a$swap_reference, length(w), na.keep, incomp.rm)

  new <- new("bed.matrix", bed = bed, snps = x@snps[w,], ped = x@ped,
           p = NULL, mu = NULL, sigma = NULL,
           standardize_p = FALSE, standardize_mu_sigma = FALSE )

  if(getOption("gaston.auto.set.stats", TRUE)) 
    new <- set.stats.snps(new, verbose = FALSE)

  w <- sum(a$swap_reference)
  if(w > 0) 
    warning(w, " reference allele inversions were performed to remove duplicated SNP")
  w <- sum(a$flip_strand)
  if(w > 0)
    warning(w, " allele strand flips were performed to remove duplicated SNP")
  if(a$NAs > 0)
  {
    if (incomp.rm) warning(a$NAs, " duplicated SNPs were removed because their alleles were incompatible")
	else warning(a$NAs, ' duplicated SNPs have incompatible alleles. Use "SNP.duplicated" function to identify them.')
  }

  new
}
