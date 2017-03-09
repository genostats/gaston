set.hwe <- function(x, method = c("chisquare", "exact"), verbose = getOption("gaston.verbose",TRUE)) {
  if( !all(c("N0", "N1", "N2") %in% names(x@snps) )) {
    if(verbose) cat("Computing basic stats\n")
    x <- set.stats(x)
  }

  w.a <- is.autosome(x@snps$chr)
  w.x <- is.chr.x(x@snps$chr)

  hwe <- rep(1.0, ncol(x) )
  method <- match.arg(method)
  if(method == 'chisquare') {
    if(verbose) cat("Computing HW chi-square p-values\n")
    hwe[w.a] <- .Call('gg_hwe_chi', PACKAGE = "gaston", x@snps$N0[w.a], x@snps$N1[w.a], x@snps$N2[w.a])
    hwe[w.x] <- .Call('gg_hwe_chi', PACKAGE = "gaston", x@snps$N0.f[w.x], x@snps$N1.f[w.x], x@snps$N2.f[w.x])
  } else {
    if(verbose) cat("Computing HW exact test p-values\n")
    hwe[w.a] <- .Call('gg_hwe', PACKAGE = "gaston", x@snps$N0[w.a], x@snps$N1[w.a], x@snps$N2[w.a])
    hwe[w.x] <- .Call('gg_hwe', PACKAGE = "gaston", x@snps$N0.f[w.x], x@snps$N1.f[w.x], x@snps$N2.f[w.x])
  }
  x@snps$hwe <- hwe
  x
}



