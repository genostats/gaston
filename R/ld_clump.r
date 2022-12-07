
LD.clump <- function(x, p, r2.threshold, p.threshold, max.dist = 500e3) {
  if(!is(x, "bed.matrix")) 
    stop("x is not a bed matrix")
  if(is.data.frame(p)) {
    if(!all(p$chr == x@ped$chr, na.rm = TRUE) | !all(p$pos == x@ped$pos, na.rm = TRUE))
      stop("Unmatching SNPs")
    a <- p
    p <- p$p
  } else {
    a <- NULL
  }
  if(length(p) != ncol(x))
    stop("Dimensions mismatch")

  or <- order(p) - 1L # -1 pour des indices qui démarrent à 0 ds la fction c++
  if(!missing(p.threshold)) {
    or <- or[ p[or + 1L] < p.threshold ]
    if(length(or) == 0) stop("too stringent p.threshold ?")
  }

  I <- .Call("ld_clump", PACKAGE = "gaston", x@bed, x@mu, x@sigma, r2.threshold, x@snps$pos, x@snps$chr, max.dist, or);

  if(is.null(a))
    a <- data.frame( chr = x@snps$chr, id = x@snps$id, pos = x@snps$pos, p = p, cluster = I+1)
  else {
    a$cluster <- I+1
  }
  a$is.index <- (a$cluster == 1:nrow(a))
  a
}


