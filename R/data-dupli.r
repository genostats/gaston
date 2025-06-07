#' Small data set to illustrate \code{\link{SNP.rm.duplicates}}
#' @name dupli
#' @aliases dupli dupli.gen dupli.ped dupli.bim dupli.pop
#' @docType data
#' 
#' @description The SNPs in this data frame are as follows:
#' \describe{
#' \item{SNP 1.}{ Unduplicated SNP }
#' \item{SNPs 2a and 2b.}{ Two duplicated SNPs with identical alleles }
#' \item{SNPs 3a and 3b.}{ Two duplicated SNPs with swapped alleles }
#' \item{SNPs 4a and 4b.}{ Two duplicated SNPs with flipped reference strand }
#' \item{SNPs 5a and 5b.}{ Two duplicated SNPs with swapped alleles and flipped reference strand }
#' \item{SNPs 6a and 6b.}{ Two duplicated SNPs with incompatible alleles }
#' \item{SNPs 7a and 7b.}{ Two duplicated SNPs including one monomorphic SNP (one allele set to \code{"0"}) }
#' \item{SNPs 8a, 8b and 8c.}{ Three duplicated SNPs }
#' \item{SNPs 9a, 9b and 9c.}{ Three duplicated SNPs with incompatible alleles }
#' }
#' 
#' @format
#' There are three data objects in the dataset:
#' \describe{
#' \item{list("dupli.gen")}{ Genotype matrix }
#' \item{list("dupli.ped")}{ Data frame containing all variables corresponding to a \code{.fam} file  }
#' \item{list("dupli.bim")}{ Data frame containing all variables corresponding to a \code{.bim} file }
#' }
#' 
#' @seealso  \code{\link{SNP.rm.duplicates}}
#' 
#' @keywords datasets
#' @examples
#' 
#' data(dupli)
#' x <- as.bed.matrix(dupli.gen, fam = dupli.ped, bim = dupli.bim)
#' 
NULL

