#' Estimating consanguinity levels from runs of homozygosity
#' 
#' @description Estimates consanguinity levels as the proportion of the genome in runs of homozygosity
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
#' @details See `\link{ROH}` for details on the definition of runs of hmozygosity
#' and for a detailed description of the parameters.
#' This functions estimates the coefficient of consanguinity \eqn{f} of each individual
#' in `x` as the proportion of its genome in runs of homozygosity.
#'
#' @return A data frame with a line for each individual in `x`, and columns `nbSNPs` 
#' (the number SNPs in ROHs), `nbSegments` (the number of ROH segments)
#' `length` (the total length of ROH segments) and `fROH` (the estimated coefficient of
#' consanguinity).
#'
#' @seealso `\link{ROH}`
#'
#' @export
fROH <- function(x, dist.units = c("bases", "cM"), NAs = c("het", "hom"), minNbSNPs = 100L, minROHlength, minDistHet, maxGapLength, beg = 1L, end = ncol(x)) {
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
    check.dist(x) # check if dist is set
    positions <- x@snps$dist
    if(missing(minROHlength)) 
      minROHlength <- 1
    if(missing(minDistHet)) 
      minDistHet <- 0.5
    if(missing(maxGapLength)) 
      maxGapLength <- 1
  }
  NAs <- match.arg(NAs)
  A <- ROHlen(x@bed, x@snps$chr, positions, beg - 1L, end - 1L, minNbSNPs, minROHlength, minDistHet, maxGapLength, NAs == "het")
  # compute total length
  CHR <- x@snps$chr[beg:end]
  chrs <- unique(CHR)
  total.len <- 0
  for(c in chrs) {
    total.len <- total.len + diff( range(positions[ CHR == c ]) )
  }
  A$fROH <- A$length / total.len
  data.frame(A)
}
