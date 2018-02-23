########################### remove duplicated

SNP.duplicated.rm <- function(x, by='chr:pos', na.keep=TRUE) {
  if(class(x) != "bed.matrix") stop("x be a bed matrix")
  
  b <- strsplit(by,':')[[1]]
  if ('alleles' %in% b) b <- c(b[b!='alleles'], 'A1', 'A2')
  
  # where are duplicated
  dupli <- SNP.match(x@snps[,b], x@snps[SNP.duplicated(x, by=by),])$index
	
  a <- .Call("gg_alleles_duplicated",  PACKAGE = "gaston", x@snps, dupli)
  bed <- .Call("gg_duplicated_remove",  PACKAGE = "gaston", x@bed, dupli, a$keep, a$flip, sum(a$keep), na.keep)

  new <- new("bed.matrix", bed = bed, snps = x@snps[a$keep,], ped = x@ped,
           p = NULL, mu = NULL, sigma = NULL,
           standardize_p = FALSE, standardize_mu_sigma = FALSE )

  if(getOption("gaston.auto.set.stats", TRUE)) 
    new <- set.stats.snps(new, verbose = FALSE)

  if(a$swap_reference > 0) 
    warning(a$swap_reference, " reference allele inversions were performed during removing duplicated SNP")
  if(a$flip_strand > 0)
    warning(a$flip_strand, " allele strand flips were performed during  removing duplicated SNP")
  if(a$NAs > 0) 
    warning(a$NAs, " duplicated SNPs were set to NA because alleles were incompatibles")

  new
}
