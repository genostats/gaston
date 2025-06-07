#' Create a \code{\link{bed.matrix}} from VCF files
#' 
#' @description  Create a \code{\link{bed.matrix}} from a \code{.vcf} file.
#' 
#' @param file  The name of the VCF file to read
#' @param max.snps  The maximal number of SNPs to read
#' @param get.info  If \code{TRUE}, the INFO field from the VCF file will integrated
#' in \code{@ped$info}
#' @param convert.chr  If \code{TRUE}, chromosomes ids \code{"X"}, \code{"Y"} and \code{"MT"}
#' will be converted in their numeric equivalents
#' @param verbose  If \code{TRUE}, display information on the function progress
#' 
#' @details
#' The vcf format is described in \url{https://github.com/samtools/hts-specs}
#' 
#' In addition to the usual data in the slot \code{@snps}, the bed.matrices produced by \code{read.vcf} have
#' \code{@snps$quality} and \code{@snps$filter} columns corresponding to the QUAL and FILTER fields in the VCF
#' file. If \code{get.info = TRUE}, an additionnal column \code{@snps$info} is added, corresponding to the
#' INFO field.
#' 
#' The information about individuals in VCF files is incomplete: in the slot \code{@ped}, the columns
#' \code{@ped$famid} and \code{@ped$id} will both contain the sample id; sex and phenotypes will be set
#' to unknown.
#' 
#' The function currently assumes that the \code{GT} field is the first field in the genotypes format.
#' If it is not the case, the variants are discarded.
#' 
#' @return  A \code{\link{bed.matrix}}
#' @seealso  \code{\link{read.bed.matrix}}
#' 
#' @keywords  vcf files
#' @examples
#' 
#' ## Read vcf file (from file name)
#' filepath <-system.file("extdata", "LCT.vcf.gz", package="gaston")
#' x1 <- read.vcf( filepath )
#' x1
#' 
#' @export read.vcf
read.vcf <- function(file, max.snps, get.info = FALSE, convert.chr = TRUE, verbose = getOption("gaston.verbose",TRUE)) {
  xx <- NULL;
  filename <- path.expand(file)

  if(missing(max.snps)) max.snps = -1L;

  L <- .Call(`_gaston_read_vcf2`, PACKAGE = "gaston", filename, max.snps, get.info)

  snp <- data.frame(chr = L$chr, id = L$id, dist = 0, pos = L$pos , A1 = L$A1, A2 = L$A2, 
                    quality = L$quality, filter = factor(L$filter), stringsAsFactors = FALSE)
  if(get.info) snp$info <-  L$info
  if(convert.chr) {
    chr <- as.integer(L$chr)
    chr[L$chr == "X"  | L$chr == "x"]  <- getOption("gaston.chr.x")[1]
    chr[L$chr == "Y"  | L$chr == "y"]  <- getOption("gaston.chr.y")[1]
    chr[L$chr == "MT" | L$chr == "mt"] <- getOption("gaston.chr.mt")[1]
    if(any(is.na(chr))) 
      warning("Some unknown chromosomes id's (try to set convert.chr = FALSE)")
    snp$chr <- chr
  } 

  ped <- data.frame(famid = L$samples, id = L$samples, father = 0, mother = 0, sex = 0, pheno = NA, stringsAsFactors = FALSE)
  x <- new("bed.matrix", bed = L$bed, snps = snp, ped = ped,
           p = NULL, mu = NULL, sigma = NULL, standardize_p = FALSE,
           standardize_mu_sigma = FALSE )

  if(getOption("gaston.auto.set.stats", TRUE)) x <- set.stats(x, verbose = verbose)
  x
}

