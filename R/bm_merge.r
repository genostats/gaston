# TODO adapter pour prendre une liste

bed.merge.inds <- function(x1, x2) {
  # on construit ped d'abord
  ped <- rbind(x1@ped, x2@ped)
  w <- !duplicated( ped[,c("famid", "id")] )
  ped <- ped[w,]
  w <- order(ped$famid, ped$id)
  ped <- ped[w,]

  x1_ids <- paste( x1@ped$famid, x1@ped$id, sep = rawToChar(as.raw(10)) ) # newline can't interfer with chars in ids!
  x2_ids <- paste( x2@ped$famid, x2@ped$id, sep = rawToChar(as.raw(10)) )
  ids    <- paste( ped$famid, ped$id, sep = rawToChar(as.raw(10)) )

  w1 <- match(ids, x1_ids, 0L)
  w2 <- match(ids, x2_ids, 0L)

  bed1 <- .Call('gg_extract_inds_indices', x1@bed, w1) 
  bed2 <- .Call('gg_extract_inds_indices', x2@bed, w2) 

  bed <- .Call("gg_bind_snps",  PACKAGE = "gaston", list(bed1, bed2))

  snps <- rbind(x1@snps, x2@snps)

  x <- new("bed.matrix", bed = bed, snps = snps, ped = ped,
           p = NULL, mu = NULL, sigma = NULL, standardize_p = FALSE,
           standardize_mu_sigma = FALSE )
  if(getOption("gaston.auto.set.stats", TRUE)) 
    x <- set.stats(x, verbose = getOption("gaston.verbose", TRUE))
  x
}



