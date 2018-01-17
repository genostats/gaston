########################### cbind


bed_cbind <- function(..., deparse.level=1) {
  L <- list(...)
  M <- lapply(L, function(x) x@bed)

  if(!all.eq( lapply(L, function(x) x@ped$famid)) | !all.eq( lapply(L, function(x) x@ped$id)) )
    stop("Individuals famids / ids are not identical, can't bind matrices")

  common_colnames <- Reduce(intersect, lapply(L, function(x) colnames(x@snps)))
  snps <- do.call(rbind, lapply(L, function(x) x@snps[common_colnames]))
  p  <- do.call(c, lapply(L, function(x) x@p) )
  mu <- do.call(c, lapply(L, function(x) x@mu))
  si <- do.call(c, lapply(L, function(x) x@sigma))
  if(length(p)  != nrow(snps)) p <- NULL
  if(length(mu) != nrow(snps)) mu <- NULL
  if(length(si) != nrow(snps)) si <- NULL
  if(anyDuplicated(snps$id))
    warning("Duplicated SNPs id's")

  if(anyDuplicated(snps[, c("chr", "pos")]))
    warning("Duplicated SNPs positions")

  bed <- .Call("gg_bind_snps",  PACKAGE = "gaston", M)
  x <- new("bed.matrix", bed = bed, snps = snps, ped = L[[1]]@ped,
           p = p, mu = mu, sigma = si, 
           standardize_p = FALSE, standardize_mu_sigma = FALSE )
  if(getOption("gaston.auto.set.stats", TRUE))
    x <- set.stats.ped(x, verbose = FALSE)
  x
}


setGeneric("cbind", signature="...")
setMethod("cbind", signature=c(...="bed.matrix"), definition = bed_cbind)
