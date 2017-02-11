read.vcf2 <- function(filename, max.snps, get.info = FALSE, convert.chr = TRUE, verbose = getOption("gaston.verbose",TRUE)) {
  filename <- path.expand(filename)
  xx <- WhopGenome::vcf_open(filename)
  if(is.null(xx)) stop("File not found")
  samples <- WhopGenome::vcf_getsamples(xx) 
  WhopGenome::vcf_selectsamples( xx, samples )

  f <- function() WhopGenome::vcf_readLineRaw(xx) 

  if(missing(max.snps)) max.snps = -1L;

  if(verbose) cat("Reading diallelic variants for", length(samples), "individuals\n");
  
  L <- .Call("gg_read_vcf2", PACKAGE="gaston", f, length(samples), max.snps, get.info)
  WhopGenome::vcf_close(xx)

  snp <- data.frame(chr = L$chr, id = L$id, dist = 0, pos = L$pos , A1 = L$A1, A2 = L$A2, 
                    quality = L$quality, filter = factor(L$filter), stringsAsFactors = FALSE)
  if(get.info) snp$info <-  L$info
  if(convert.chr) {
    chr <- as.numeric(L$chr)
    chr <- ifelse(L$chr == "X"  | L$chr == "x",  options("gaston.chr.x")[1],  chr)
    chr <- ifelse(L$chr == "Y"  | L$chr == "y",  options("gaston.chr.y")[1],  chr)
    chr <- ifelse(L$chr == "MT" | L$chr == "mt", options("gaston.chr.mt")[1], chr)
    if(any(is.na(chr))) 
      warning("Some unknown chromosomes id's (try to set convert.chr = FALSE)")
    snp$chr <- chr
  } 

  ped <- data.frame(famid = samples, id = samples, father = 0, mother = 0, sex = 0, pheno = NA, stringsAsFactors = FALSE)
  x <- new("bed.matrix", bed = L$bed, snps = snp, ped = ped,
           p = NULL, mu = NULL, sigma = NULL, standardize_p = FALSE,
           standardize_mu_sigma = FALSE )

  if(getOption("gaston.auto.set.stats", TRUE)) x <- set.stats(x, verbose = verbose)
  x
}


read.vcf <- function(filename, max.snps, verbose = getOption("gaston.verbose",TRUE)) {
  filename <- path.expand(filename)
  xx <- WhopGenome::vcf_open(filename)
  if(is.null(xx)) stop("File not found")
  samples <- WhopGenome::vcf_getsamples(xx) 
  WhopGenome::vcf_selectsamples( xx, samples )

  f <- function() WhopGenome::vcf_readLineRaw(xx) 

  if(missing(max.snps)) {
    if(verbose) cat("Counting diallelic variants\n")
    max.snps <- .Call("gg_count_dia_vcf", PACKAGE="gaston", f)
    WhopGenome::vcf_close(xx);
    WhopGenome::vcf_reopen(xx);
  }

  if(verbose) cat("Reading", max.snps, "diallelic variants for", length(samples), "individuals\n");
  
  L <- .Call("gg_read_vcf", PACKAGE="gaston", f, length(samples), max.snps)
  WhopGenome::vcf_close(xx)
  snp <- data.frame(chr = L$chr, id = L$id, dist = 0, pos = L$pos , A1 = L$A1, A2 = L$A2, stringsAsFactors = FALSE)
  ped <- data.frame(famid = samples, id = samples, father = 0, mother = 0, sex = 0, pheno = NA, stringsAsFactors = FALSE)
  x <- new("bed.matrix", bed = L$bed, snps = snp, ped = ped,
           p = NULL, mu = NULL, sigma = NULL, standardize_p = FALSE,
           standardize_mu_sigma = FALSE )

  if(getOption("gaston.auto.set.stats", TRUE)) x <- set.stats(x, verbose = verbose)
  x
}

count.vcf <- function(filename) {
  xx <- WhopGenome::vcf_open(filename)
  if(is.null(xx)) stop("File not found")
  samples <- WhopGenome::vcf_getsamples(xx) 
  WhopGenome::vcf_selectsamples( xx, samples )
  f <- function() WhopGenome::vcf_readLineRaw(xx) 
  n <- .Call("gg_count_dia_vcf", PACKAGE="gaston", f)
  WhopGenome::vcf_close(xx)
  n
}


