# NAs = should missing data be treated as heterozygous or homozygous ?
# minROHlength = minimal length of a ROH
# minDistHet   = minimal distance between hetorozygous SNPs
# maxGapLength = maximal distance between two SNPs in a ROH

fROH <- function(x, dist.units = c("bases", "cM"), NAs = c("het", "hom"), minROHlength, minDistHet, maxGapLength, beg = 1L, end = ncol(x)) {
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
    positions <- x@snps$dist
    if(missing(minROHlength)) 
      minROHlength <- 1
    if(missing(minDistHet)) 
      minDistHet <- 0.5
    if(missing(maxGapLength)) 
      maxGapLength <- 1
  }
  NAs <- match.arg(NAs)
  A <- ROHlen(x@bed, x@snps$chr, positions, beg - 1L, end - 1L, minROHlength, minDistHet, maxGapLength, NAs == "het")
  # compute total length
  chrs <- unique(x@snps$chr)
  total.len <- 0
  for(c in chrs) {
    total.len <- total.len + diff( range(positions[ x@snps$chr == c ]) )
  }
  A$fROH <- A$length / total.len
  data.frame(A)
}
