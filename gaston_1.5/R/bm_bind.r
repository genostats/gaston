########################### cbind

setGeneric("cbind", signature="...")
setMethod("cbind", signature=c(...="bed.matrix"), definition= function(..., deparse.level=1) {
  L <- list(...)
  M <- lapply(L, function(x) x@bed)

  if(!all.eq( lapply(L, function(x) x@ped$famid)) | !all.eq( lapply(L, function(x) x@ped$id)) )
    stop("Individuals famids / ids are not identical, can't bind matrices")

  common_colnames <- Reduce(intersect, lapply(L, function(x) colnames(x@snps)))
  snps <- Reduce(rbind, lapply(L, function(x) x@snps[common_colnames]))

  if(anyDuplicated(snps$id))
    warning("Duplicated SNPs id's")

  if(anyDuplicated(snps[, c("chr", "pos")]))
    warning("Duplicated SNPs positions")

  bed <- .Call("gg_bind_snps",  PACKAGE='gaston', M)
  x <- new("bed.matrix", bed = bed, snps = snps, ped = L[[1]]@ped,
           p = NULL, mu = NULL, sigma = NULL, 
           standardize_p = FALSE, standardize_mu_sigma = FALSE )
  if(getOption("gaston.auto.set.stats", TRUE))
    x <- set.stats.ped(x, verbose = FALSE)
  x
})









########################### rbind

setGeneric("rbind", signature="...")
setMethod("rbind", signature=c(...="bed.matrix"), definition= function(..., deparse.level=1) {
  L <- list(...)
  M <- lapply(L, function(x) x@snps)

  if(!all.eq( lapply(M, function(x) x$id))) 
    stop("SNP ids are not identical, can't bind matrices")

  common_colnames <- Reduce(intersect, lapply(L, function(x) colnames(x@ped)))
  ped <- Reduce(rbind, lapply(L, function(x) x@ped[common_colnames]))

  if(anyDuplicated(ped[, c("famid", "id")]))
    warning("There are duplicated individuals (same family and individual id)")

  a <- .Call("gg_alleles_recoding",  PACKAGE='gaston', M)
  M <- lapply(L, function(x) x@bed)
  bed <- .Call("gg_bind_inds2",  PACKAGE='gaston', M, a$flip)

  x <- new("bed.matrix", bed = bed, snps = L[[1]]@snps, ped = ped,
           p = NULL, mu = NULL, sigma = NULL,
           standardize_p = FALSE, standardize_mu_sigma = FALSE )

  if(getOption("gaston.auto.set.stats", TRUE)) 
    x <- set.stats.snps(x, verbose = FALSE)

  if(a$ref > 0) 
    warning(a$ref, " reference allele inversions were performed during binding")
  if(a$strand > 0)
    warning(a$strand, " allele strand flips were performed during binding")
  if(a$NAs > 0) 
    warning(a$NAs, " SNPs were set to NA because alleles were incompatibles")

  x
})

# pour vérifier les rs id... / famid id
# le as.character() sert à éviter des soucis si il y a des facteurs
all.eq <- function(L) {
  a <- lapply(L[-1], function(x) all(as.character(x) == as.character(L[[1]])))
  Reduce("&" ,a)
}

# keep that private, for testing only
alleles.recoding <- function(...) {
  L <- list(...)
  M <- lapply(L, function(x) x@snps)
  .Call("gg_alleles_recoding",  PACKAGE='gaston', M)
}


# ! inplace modifications
invert_snp_coding <- function(x, snp) {
  .Call("gg_invert_snp_coding",  PACKAGE='gaston', x@bed, snp)
}

snp_hz_to_na <- function(x, snp) {
  .Call("gg_snp_hz_to_na",  PACKAGE='gaston', x@bed, snp)
}

set_snp_to_na <- function(x, snp) {
  .Call("gg_set_snp_to_na",  PACKAGE='gaston', x@bed, snp)
}

