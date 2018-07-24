########################### remove duplicated

SNP.rm.duplicates <- function(x, by = "chr:pos", na.keep = TRUE, incomp.rm = TRUE) {
  if(class(x) != "bed.matrix") stop("x should be a bed.matrix")
  
  b <- strsplit(by,':')[[1]]
  if ('alleles' %in% b) b <- c(b[b!='alleles'], 'A1', 'A2')
  
  # where are duplicated SNPs
  dupli <- SNP.match(x@snps[,b,drop=FALSE], x@snps[SNP.duplicated(x, by=by),], by = by)$index
	
  a <- .Call("gg_alleles_duplicated",  PACKAGE = "gaston", x@snps, dupli)
  
  if (incomp.rm) w <- which(a$keep == TRUE)
  else  w <- which(is.na(a$keep) | a$keep == TRUE)
  
  bed <- .Call("gg_duplicated_remove",  PACKAGE = "gaston", x@bed, dupli, a$keep, a$swap_reference, length(w), na.keep, incomp.rm)

  new <- new("bed.matrix", bed = bed, snps = x@snps[w,], ped = x@ped,
           p = NULL, mu = NULL, sigma = NULL,
           standardize_p = FALSE, standardize_mu_sigma = FALSE )

  if(getOption("gaston.auto.set.stats", TRUE)) 
    new <- set.stats.snps(new, verbose = FALSE)

  w <- sum(a$swap_reference)
  if(w > 0) 
    warning(w, " reference allele inversions were performed to remove duplicated SNP")
  w <- sum(a$flip_strand)
  if(w > 0)
    warning(w, " allele strand flips were performed to remove duplicated SNP")
  if(a$NAs > 0)
  {
    if (incomp.rm) warning(a$NAs, " duplicated SNPs were removed because their alleles were incompatible")
	else warning(a$NAs, ' duplicated SNPs have incompatible alleles. Use "SNP.duplicated" function to identify them.')
  }

  new
}
