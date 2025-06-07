#' @rdname set.stats
#' @export
set.stats.ped <- function(x, verbose = getOption("gaston.verbose",TRUE)) {
  if (!is(x, "bed.matrix")) stop('x must be an object of class bed.matrix')

  w.a  <- is.autosome(x@snps$chr)
  w.x  <- is.chr.x(x@snps$chr)
  w.y  <- is.chr.y(x@snps$chr)
  w.mt <- is.chr.mt(x@snps$chr)

  st <- .Call(`_gaston_geno_stats_inds`, PACKAGE = "gaston", x@bed, w.x, w.y, w.mt) 

  ############ completer inds/ped
  n.a <- sum(w.a)
  st$inds$callrate <- 1-st$inds$NAs/n.a
  st$inds$hz <- st$inds$N1/(n.a-st$inds$NAs)

  n.x <- sum(w.x)
  st$inds$callrate.x <- 1-st$inds$NAs.x/n.x
  st$inds$hz.x <- st$inds$N1.x/(n.x-st$inds$NAs.x)

  n.y <- sum(w.y)
  st$inds$callrate.y <- 1-st$inds$NAs.y/n.y
  st$inds$hz.y <- st$inds$N1.y/(n.y-st$inds$NAs.y)

  n.mt <- sum(w.mt)
  st$inds$callrate.mt <- 1-st$inds$NAs.mt/n.mt
  st$inds$hz.mt <- st$inds$N1.mt/(n.mt-st$inds$NAs.mt)

  ########## insérer dans x
  x@ped[ , names(st$inds)] <- st$inds

  if(verbose) cat("ped stats have been set. \n") 
  x
}

