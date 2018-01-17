########################### rbind

bed_rbind <- function(..., deparse.level = 1) {
  L <- list(...)
  M <- lapply(L, function(x) x@snps)

  # ici il faudrait faire un truc un peu plus fin ... (plutôt vérifier chr:pos:alleles...)
  if(!all.eq( lapply(M, function(x) x$id))) 
    stop("SNP ids are not identical, can't bind matrices")

  common_colnames <- Reduce(intersect, lapply(L, function(x) colnames(x@ped)))
  ped <- Reduce(rbind, lapply(L, function(x) x@ped[common_colnames]))

  if(anyDuplicated(ped[, c("famid", "id")]))
    warning("There are duplicated individuals (same family and individual id)")

  a <- .Call("gg_alleles_recoding",  PACKAGE = "gaston", M)
  M <- lapply(L, function(x) x@bed)
  bed <- .Call("gg_bind_inds2",  PACKAGE = "gaston", M, a$flip)

  x <- new("bed.matrix", bed = bed, snps = L[[1]]@snps, ped = ped,
           p = NULL, mu = NULL, sigma = NULL,
           standardize_p = FALSE, standardize_mu_sigma = FALSE )

  if(getOption("gaston.auto.set.stats", TRUE)) 
    x <- set.stats.snps(x, verbose = FALSE)

  if(a$swap_reference > 0) 
    warning(a$swap_reference, " reference allele inversions were performed during binding")
  if(a$flip_strand > 0)
    warning(a$flip_strand, " allele strand flips were performed during binding")
  if(a$NAs > 0) 
    warning(a$NAs, " SNPs were set to NA because alleles were incompatibles")

  x
}

setGeneric("rbind", signature="...")
setMethod("rbind", signature=c(...="bed.matrix"), definition = bed_rbind)
