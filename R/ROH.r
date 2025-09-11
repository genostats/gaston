#' Runs of Homozygosity
#' 
#' @description Identifies Runs of Homozygosity
#'
#' @param x a bed.matrix
#' @param dist.units distance unit to use (either 'bases' or 'cM')
#' @param NAs Should missing data be treated as heterozygous or homozygous?
#' @param minNbSNPs minimal number of SNPS in a ROH
#' @param minROHlength minimal length of a ROH (in dist.units)
#' @param minDistHet minimal distance between hetorozygous SNPs (in dist.units)
#' @param maxGapLength maximal distance between two SNPs in a ROH (in dist.units=
#' @param beg index of the first SNP to consider
#' @param end index of the last SNP to consider
#'
#' @details Runs of homozygosity (ROHs) are sequences of homozygous SNPs, including a minimal number 
#' `minNbSNPs` of SNPs, and spanning a minimal length of `minROHlength`, in units specified
#' by `dist.units`. They can contain some heterozygous SNPs, if these are at distance 
#' larger than `minDistHet` from the extremities of the ROH, and from each other.
#' Consecutive SNPs in a ROH are at a distance smaller than `maxGepLength`.
#'
#' The distance units are either base pairs ("bases") or centiMorgans ("cM").
#' It is better to use `dist.units = "cM"` but this requires that the `dist`
#' column in `x@snps` is set (see `\link{set.dist}`).
#'
#' If `minROHlength`, `minDistHet` and `maxGapLength` are not specified, default values
#' are set depending on the units chosen as follows: if `dist.units = "bases"` 
#' `minROHlength = 1e6` (1Mb), `minDistHet = 5e5` (500 kb) and `maxGapLength = 1e6`.
#' If `dist.units = "cM"`, `minROHlength = 1`, `minDistHet = 0.5` and `maxGapLength = 1`.
#'
#' @return A list of data frames, one for each indivual in `x`, describing the ROHs.
#' These data frames have columns `i.beg`, the index of the first SNPs in the ROH,
#' `i.end`, the index of the last SNP in the ROH, `pos.beg` and `pos.end` their genomic
#' position (in accordance with `dist.units`) and `nb.het`, the number of heterozygous
#' SNPs in the ROH.
#'
#' @seealso `\link{set.dist}`, `\link{fROH}`
#'
#' @export
ROH <- function(x, dist.units = c("bases", "cM"), NAs = c("het", "hom"), minNbSNPs = 100L, minROHlength, minDistHet, maxGapLength, beg = 1L, end = ncol(x)) {
  dist.units <- match.arg(dist.units)
  if(dist.units == "bases") {
    positions <- x@snps$pos
    if(missing(minROHlength)) 
      minROHlength <- 1e6
    if(missing(minDistHet)) 
      minDistHet <- 0.5e6
    if(missing(maxGapLength)) 
      maxGapLength <- 1e6
  } else {
    check.dist(x) # check if 'dist' is set
    positions <- x@snps$dist
    if(missing(minROHlength)) 
      minROHlength <- 1
    if(missing(minDistHet)) 
      minDistHet <- 0.5
    if(missing(maxGapLength)) 
      maxGapLength <- 1
  }
  NAs <- match.arg(NAs)
  ROH_(x@bed, x@snps$chr, positions, beg - 1L, end - 1L, minNbSNPs, minROHlength, minDistHet, maxGapLength, NAs == "het")
}
