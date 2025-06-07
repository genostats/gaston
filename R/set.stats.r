#' Basic statistics for a \code{\link{bed.matrix}}
#' 
#' @description
#' Return an updated \code{\link{bed.matrix}} with new variables for
#' several basic statistics.
#' 
#' @param x  A \code{\link{bed.matrix}}
#' @param set.p  If \code{TRUE}, \code{x@p} is updated
#' @param set.mu_sigma  If \code{TRUE}, \code{x@mu} and \code{x@sigma} are updated
#' @param verbose  If \code{TRUE}, display information on the function actions
#' 
#' 
#' @details
#' \code{set.stats} is called by default by all functions that create a bed.matrix, unless
#' the global option \code{gaston.auto.set.stats} is \code{FALSE} (cf example below).
#' 
#' \code{set.stats} and \code{set.stats.ped} update \code{x@ped}, adding the following variables:
#' \itemize{
#' \item \code{N0}, \code{N1}, \code{N2} and \code{NAs} give for each individual the number of
#' autosomal SNPs with a genotype equal to 0, 1, 2 and missing, respectively
#' \item \code{N0.x}, \code{N1.x}, \code{N2.x} and \code{NAs.x} idem for chromosome X
#' \item \code{N0.y}, \code{N1.y}, \code{N2.y} and \code{NAs.y} idem for chromosome Y
#' \item \code{N0.mt}, \code{N1.mt}, \code{N2.mt} and \code{NAs.mt} idem for mitochondrial SNPs
#' \item \code{callrate}, \code{callrate.x}, \code{callrate.y}, \code{callrate.mt}
#' is the individual callrate for autosomal, X, Y, mitochondrial SNPs
#' \item \code{hz}, \code{hz.x}, \code{hz.y}, \code{hz.mt} is the individual heterozygosity
#' for autosomal, X, Y, mitochondrial SNPs
#' }
#' 
#' \code{set.stats} and \code{set.stats.snps} update \code{x@snps}, adding the following variables:
#' \itemize{
#' \item \code{N0}, \code{N1}, \code{N2} and \code{NAs} give for each SNP the number of individuals
#' with a genotype equal to 0, 1, 2 and missing, respectively
#' \item \code{N0.f}, \code{N1.f}, \code{N2.f} and \code{NAs.f} give, only for SNPs on chromosome X,
#' the number of female individuals with a genotype equal to 0, 1, 2 and missing, respectively
#' \item \code{callrate} is the SNP callrate (for Y linked SNPs, the callrate is computed usin
#' males only).
#' \item \code{maf} is the Minor Allele Frequency
#' \item \code{hz} is the SNP heterozygosity (for X linked SNPs, the heterozygosity is computed
#' using females only).
#' }
#' 
#' If \code{set.p = TRUE}, \code{x@p} is updated with the alternate allele frequency.
#' 
#' If \code{set.mu_sigma = TRUE}, \code{x@mu} is updated with the genotype mean (equal to \code{2*x@p})
#' and \code{x@sigma} with the genotype standard deviation (should be approximately \code{sqrt(2*x@p*(1-x@p))}
#' under Hardy-Weinberg Equilibrium).
#' 
#' @return
#' A \code{\link{bed.matrix}} similar to \code{x}, with slots updated as described above.
#' @seealso  \code{\link{set.hwe}}, \code{\link{set.genomic.sex}}
#' 
#' @keywords  Genetic   Statistics
#' @examples
#' 
#' # Disable auto set stats :
#' options(gaston.auto.set.stats = FALSE)
#' 
#' # Load data
#' data(TTN)
#' x <- as.bed.matrix(TTN.gen, TTN.fam, TTN.bim)
#' str(x@ped)
#' str(x@snps)
#' 
#' # Compute statistics
#' x <- set.stats(x)
#' str(x@ped)
#' str(x@snps)
#' 
#' # restore default behavior
#' options(gaston.auto.set.stats = TRUE)
#' 
#' @export set.stats
set.stats <- function(x, set.p = TRUE, set.mu_sigma = TRUE, verbose = getOption("gaston.verbose",TRUE)) {
  if (!is(x, "bed.matrix")) stop('x must be an object of class bed.matrix')
  if(!is.logical(set.p) | !is.logical(set.mu_sigma)) 
    stop('set.* arguments must be logical')

  w.a  <- is.autosome(x@snps$chr)
  w.x  <- is.chr.x(x@snps$chr)
  w.y  <- is.chr.y(x@snps$chr)
  w.mt <- is.chr.mt(x@snps$chr)
  w.f  <- x@ped$sex  == 2
  w.f[is.na(w.f)] <- FALSE   # tout ce qui n'est pas sex == 2 est pris comme un homme

  st <- .Call(`_gaston_geno_stats`, PACKAGE = "gaston", x@bed, w.x, w.y, w.mt, w.f) 

  nb.f <- sum(x@ped$sex == 2)
  nb.h <- nrow(x) - nb.f
  ############  completer snps
  st$snps$callrate <- 1-st$snps$NAs/nrow(x)

  # correction pour chr y 
  # (ped$sed peut être faux / on compte les NAs chez les individus avec sex != 2...)
  st$snps$callrate[w.y] <- 1-(st$snps$NAs[w.y]-st$snps$NAs.f[w.y])/nb.h

  # freq allele alt
  n <- nrow(x) - st$snps$NAs;
  pp <- (2*st$snps$N2 + st$snps$N1)/(2*n);

  # mise à NA des stats .f ailleurs que sur x ou y
  st$snps$N0.f[!w.x & !w.y] <- NA;  st$snps$N1.f[!w.x & !w.y] <- NA
  st$snps$N2.f[!w.x & !w.y] <- NA; st$snps$NAs.f[!w.x & !w.y] <- NA

  # correction pour chr x
  a <- st$snps$N2.f[w.x] + st$snps$N1.f[w.x] + st$snps$N2[w.x]
  b <- st$snps$N0.f[w.x] + st$snps$N1.f[w.x] + st$snps$N0[w.x]
  pp[w.x] <- a/(a+b)

  st$snps$maf <- pmin(pp,1-pp)
  st$snps$hz <- st$snps$N1/n

  # correction pour chr x
  st$snps$hz[w.x] <- st$snps$N1.f[w.x]/(nb.f - st$snps$NAs.f[w.x])

  ############ completer inds/ped
  n.a <- sum(w.a)
  st$inds$callrate <- 1-st$inds$NAs/n.a
  st$inds$hz <- st$inds$N1/(n.a-st$inds$NAs)

  n.x <- sum(w.x)
  st$inds$callrate.x <- 1-st$inds$NAs.x/n.x
  st$inds$hz.x <- st$inds$N1.x/(n.x-st$inds$NAs.x)

  n.y <- sum(w.y)
  st$inds$callrate.y <- 1-st$inds$NAs.y/n.y
  st$inds$hz.y <- st$inds$N1.y/(n.y-st$inds$NAs.y)

  n.mt <- sum(w.mt)
  st$inds$callrate.mt <- 1-st$inds$NAs.mt/n.mt
  st$inds$hz.mt <- st$inds$N1.mt/(n.mt-st$inds$NAs.mt)

  ########## insérer dans x
  x@snps[, names(st$snps)] <- st$snps
  x@ped[ , names(st$inds)] <- st$inds

  if(verbose) cat("ped stats and snps stats have been set. \n") 

  if(set.p) {
    x@p <- pp
    if(verbose) cat("'p' has been set. \n")
  }  

  if(set.mu_sigma) { # calcul brutal
    n <- nrow(x) - x@snps$NAs;
    mu <- (2*x@snps$N2 + x@snps$N1)/n; # c'est 2 pp ... enfin bref
    N <- nrow(x)
    s <- sqrt( (x@snps$N1 + 4*x@snps$N2 + mu**2*x@snps$NAs)/(N-1) - N/(N-1)*mu**2 )
    x@mu <- mu;
    x@sigma <- s
    if(verbose) cat("'mu' and 'sigma' have been set.\n");
  }
  x
}

