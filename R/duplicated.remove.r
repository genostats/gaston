########################### remove duplicated

SNP.duplicated.rm <- function(x, by='chr:pos') {
  if(class(x) != "bed.matrix") stop("x be a bed matrix")

  # where are duplicated
  dupli <- SNP.match(x@snps[strsplit(by,':')[[1]]], x@snps[SNP.duplicated(x, by=by),])
	
  a <- .Call("gg_alleles_duplicated",  PACKAGE = "gaston", x@snps, dupli)
  bed <- .Call("gg_duplicated",  PACKAGE = "gaston", bed, dupli, a$flip, sum(a$flip=='remove'))

  x <- new("bed.matrix", bed = bed, snps = bed@snps[,a$flip=='keep'], ped = x@ped$ped,
           p = NULL, mu = NULL, sigma = NULL,
           standardize_p = FALSE, standardize_mu_sigma = FALSE )

  if(getOption("gaston.auto.set.stats", TRUE)) 
    x <- set.stats.snps(x, verbose = FALSE)

  if(a$swap_reference > 0) 
    warning(a$swap_reference, " reference allele inversions were performed during removing duplicated SNP")
  if(a$flip_strand > 0)
    warning(a$flip_strand, " allele strand flips were performed during  removing duplicated SNP")
  if(a$NAs > 0) 
    warning(a$NAs, " duplicated SNPs were set to NA because alleles were incompatibles")

  x
}
