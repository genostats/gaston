# pour vérifier les rs id... / famid id
# le as.character() sert à éviter des soucis si il y a des facteurs
all_eq <- function(L) {
  a <- lapply(L[-1], function(x) all(as.character(x) == as.character(L[[1]])))
  Reduce("&" ,a)
}

# keep that private, for testing only
alleles.recoding <- function(...) {
  L <- list(...)
  M <- lapply(L, function(x) x@snps)
  .Call(`_gaston_alleles_recoding`, PACKAGE = "gaston", M)
}


# ! inplace modifications
invert_snp_coding <- function(x, snp) {
  .Call(`_gaston_invert_snp_coding`, PACKAGE = "gaston", x@bed, snp)
}

snp_hz_to_na <- function(x, snp) {
  .Call(`_gaston_snp_hz_to_na`, PACKAGE = "gaston", x@bed, snp)
}

set_snp_to_na <- function(x, snp) {
  .Call(`_gaston_set_snp_to_na`, PACKAGE = "gaston", x@bed, snp)
}

