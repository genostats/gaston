LD.thin <- function(x, threshold, max.dist = 250e3, beg = 1, end = ncol(x), which.snps,
                    dist.unit = c("bases", "indices", "cM"), extract = TRUE, keep = c("left", "right", "random")) {

  if(missing(which.snps)) which.snps <- rep(TRUE, end-beg+1)

  if(!is.logical(which.snps) | length(which.snps) != end-beg+1)
    stop("which.snps must be a Logical vector of length end - beg + 1")

  if(is.null(x@mu) | is.null(x@sigma))
    stop("LD.thin needs mu and sigma to be set for LD computation (use set.stats)")

  # ne pas considérer les SNPs monomorphes ou qui ont un callrate nul
  which.snps <- which.snps & (x@snps$callrate > 0) & (x@snps$maf > 0)

  dist.unit <- match.arg(dist.unit)
  if(dist.unit == "bases") {
    pos = as.integer(x@snps$pos)
  } else if(dist.unit == "indices") {
    pos = seq_len(ncol(x))
  } else { # cM = conversion en micro morgans !
    if(max.dist > 100)
      warning("max.dist value seems very high for dist.unit = \"cM\"")
    pos = as.integer( x@snps$dist * 1e4 )
    max.dist = as.integer(max.dist*1e4)
  }
  if( all(x@snps$pos == x@snps$pos[1]) )
    stop("Position of SNPs must be available")

  keep <- match.arg(keep)
  if(keep == "left") {
    w <- .Call("gg_ld_thin_left", x@bed, x@mu, x@sigma, threshold, pos, 
          as.integer(x@snps$chr), as.integer(max.dist), as.integer(beg)-1L, as.integer(end)-1L,
          which.snps)
  } else if (keep == "right"){
    w <- .Call("gg_ld_thin_right", x@bed, x@mu, x@sigma, threshold, pos, 
          as.integer(x@snps$chr), as.integer(max.dist), as.integer(beg)-1L, as.integer(end)-1L,
          which.snps)
  } else {
    w <- .Call("gg_ld_thin_random", x@bed, x@mu, x@sigma, threshold, pos, 
          as.integer(x@snps$chr), as.integer(max.dist), as.integer(beg)-1L, as.integer(end)-1L,
          which.snps)
  }

  if(!extract) return(w)
  x[ , seq(beg,end)[w] ]
}

