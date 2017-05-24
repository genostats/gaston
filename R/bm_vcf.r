read.vcf <- function(file, max.snps, get.info = FALSE, convert.chr = TRUE, verbose = getOption("gaston.verbose",TRUE)) {
  xx <- NULL;
  filename <- path.expand(file)

  if(missing(max.snps)) max.snps = -1L;

  L <- .Call("gg_read_vcf2", PACKAGE = "gaston", filename, max.snps, get.info)

  snp <- data.frame(chr = L$chr, id = L$id, dist = 0, pos = L$pos , A1 = L$A1, A2 = L$A2, 
                    quality = L$quality, filter = factor(L$filter), stringsAsFactors = FALSE)
  if(get.info) snp$info <-  L$info
  if(convert.chr) {
    chr <- as.integer(L$chr)
    chr[L$chr == "X"  | L$chr == "x"]  <- getOption("gaston.chr.x")[1]
    chr[L$chr == "Y"  | L$chr == "y"]  <- getOption("gaston.chr.y")[1]
    chr[L$chr == "MT" | L$chr == "mt"] <- getOption("gaston.chr.mt")[1]
    if(any(is.na(chr))) 
      warning("Some unknown chromosomes id's (try to set convert.chr = FALSE)")
    snp$chr <- chr
  } 

  ped <- data.frame(famid = L$samples, id = L$samples, father = 0, mother = 0, sex = 0, pheno = NA, stringsAsFactors = FALSE)
  x <- new("bed.matrix", bed = L$bed, snps = snp, ped = ped,
           p = NULL, mu = NULL, sigma = NULL, standardize_p = FALSE,
           standardize_mu_sigma = FALSE )

  if(getOption("gaston.auto.set.stats", TRUE)) x <- set.stats(x, verbose = verbose)
  x
}

read.vcf.filtered <- function(file, positions, max.snps, get.info = FALSE, convert.chr = TRUE, verbose = getOption("gaston.verbose",TRUE)) {
  xx <- NULL;
 filename <- path.expand(file)

  if(missing(max.snps)) max.snps = -1L;

  L <- .Call("gg_read_vcf_filtered", PACKAGE = "gaston", filename, positions, max.snps, get.info)

  snp <- data.frame(chr = L$chr, id = L$id, dist = 0, pos = L$pos , A1 = L$A1, A2 = L$A2, 
                    quality = L$quality, filter = factor(L$filter), stringsAsFactors = FALSE)
  if(get.info) snp$info <-  L$info
  if(convert.chr) {
    chr <- as.integer(L$chr)
    chr[L$chr == "X"  | L$chr == "x"]  <- getOption("gaston.chr.x")[1]
    chr[L$chr == "Y"  | L$chr == "y"]  <- getOption("gaston.chr.y")[1]
    chr[L$chr == "MT" | L$chr == "mt"] <- getOption("gaston.chr.mt")[1]
    if(any(is.na(chr))) 
      warning("Some unknown chromosomes id's (try to set convert.chr = FALSE)")
    snp$chr <- chr
  } 

  ped <- data.frame(famid = L$samples, id = L$samples, father = 0, mother = 0, sex = 0, pheno = NA, stringsAsFactors = FALSE)
  x <- new("bed.matrix", bed = L$bed, snps = snp, ped = ped,
           p = NULL, mu = NULL, sigma = NULL, standardize_p = FALSE,
           standardize_mu_sigma = FALSE )

  if(getOption("gaston.auto.set.stats", TRUE)) x <- set.stats(x, verbose = verbose)
  x
}

