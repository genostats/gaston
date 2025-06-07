#' @rdname set.stats
#' @export
set.stats.snps <- function(x, set.p = TRUE, set.mu_sigma = TRUE, verbose = getOption("gaston.verbose",TRUE)) {
  if( is(x)!='bed.matrix' ) stop('x must be an object of class bed.matrix')
  if(!is.logical(set.p) | !is.logical(set.mu_sigma)) 
    stop('set.* arguments must be logical')

  w.x  <- is.chr.x(x@snps$chr)
  w.y  <- is.chr.y(x@snps$chr)
  w.f <- x@ped$sex  == 2
  w.f[is.na(w.f)] <- FALSE   # tout ce qui n'est pas sex == 2 est pris comme un homme

  st <- .Call(`_gaston_geno_stats_snps`, PACKAGE = "gaston", x@bed, w.x|w.y, w.f) 

  nb.f <- sum(x@ped$sex == 2)
  nb.h <- nrow(x) - nb.f

  ############  completer snps
  st$snps$callrate <- 1-st$snps$NAs/nrow(x)

  # correction pour chr y 
  # (ped$sed peut être faux / on compte les NAs chez les individus avec sex != 2...)
  st$snps$callrate[w.y] <- 1-(st$snps$NAs[w.y]-st$snps$NAs.f[w.y])/nb.h


  # mise à NA des stats .f ailleurs que sur le x
  st$snps$N0.f[!w.x] <- NA; st$snps$N1.f[!w.x] <- NA
  st$snps$N2.f[!w.x] <- NA; st$snps$NAs.f[!w.x] <- NA

  # freq allele alt
  n <- nrow(x) - st$snps$NAs;
  pp <- (2*st$snps$N2 + st$snps$N1)/(2*n);

  # correction pour chr x
  a <- st$snps$N2.f[w.x] + st$snps$N1.f[w.x] + st$snps$N2[w.x]
  b <- st$snps$N0.f[w.x] + st$snps$N1.f[w.x] + st$snps$N0[w.x]
  pp[w.x] <- a/(a+b)

  st$snps$maf <- pmin(pp,1-pp)
  st$snps$hz <- st$snps$N1/n

  # correction pour chr x
  st$snps$hz[w.x] <- st$snps$N1.f[w.x]/(nb.f - st$snps$NAs.f[w.x])

  ########## insérer dans x
  x@snps[, names(st$snps)] <- st$snps

  if(verbose) cat("snps stats have been set. \n") 

  if(set.p) {
    x@p <- pp
    if(verbose) cat("'p' has been set. \n")
  }  

  if(set.mu_sigma) { # calcul brutal
    n <- nrow(x) - x@snps$NAs;
    mu <- (2*x@snps$N2 + x@snps$N1)/n;
    N <- nrow(x)
    s <- sqrt( (x@snps$N1 + 4*x@snps$N2 + mu**2*x@snps$NAs)/(N-1) - N/(N-1)*mu**2 )
    x@mu <- mu;
    x@sigma <- s
    if(verbose) cat("'mu' and 'sigma' have been set.\n");
  }
  x
}

