#' Read a \code{\link{bed.matrix}}
#' 
#' @description
#' Create a \code{\link{bed.matrix}} from a \code{.bed} file, and either
#' a \code{.rds} file or a \code{.bim} and a \code{.fam} file.
#' 
#' @param basename  Basename of all files
#' @param bed  Name of the \code{.bed} file
#' @param fam  Name of the \code{.fam} file
#' @param bim  Name of the \code{.bim} file
#' @param rds  Name of the \code{.rds} file (ignored if \code{NULL})
#' @param verbose  If \code{TRUE}, display information on the function actions
#' 
#' 
#' @details
#' The \code{.bed}, \code{.fam} and \code{.bim} files follow the PLINK specifications
#' (\url{http://zzz.bwh.harvard.edu/plink/binary.shtml}).
#' 
#' If a \code{.rds} file exists (created by \code{write.bed.matrix}),
#' the \code{.fam} and \code{.bim} files will be ignored.
#' To ignore an existing \code{.rds} file, set \code{rds = NULL}.
#' 
#' If the \code{.bed} file does not exist, and \code{basename} ends by \code{".bed"},
#' the function will try to generate a new basename by trimming the extension out. This
#' allows to write \code{read.bed.matrix("file.bed")} instead of \code{read.bed.matrix("file")}.
#' 
#' If the option \code{gaston.auto.set.stats} is set to \code{TRUE} (the default),
#' the function \code{set.stats} will be called before returning the \code{bed.matrix},
#' unless a \code{.rds} file is present: in this case, the \code{bed.matrix} obtained
#' is identical to the \code{bed.matrix} saved with \code{write.bed.matrix}.
#' 
#' 
#' @return  A \code{\link{bed.matrix}}
#' @seealso  \code{\link{write.bed.matrix}}, \code{\link{set.stats}}
#' 
#' 
#' @examples
#' 
#' # Read RDS and bed files
#' x <- read.bed.matrix( system.file("extdata", "LCT.bed", package="gaston") )
#' x
#' 
#' @export read.bed.matrix
read.bed.matrix <- function(basename, bed = paste(basename, ".bed", sep=""), fam = paste(basename, ".fam", sep=""), 
                                      bim = paste(basename, ".bim", sep=""), rds = paste(basename, ".rds", sep=""),
                            verbose = getOption("gaston.verbose",TRUE)) {

  bed <- path.expand(bed)
  if(!file.exists(bed)) { # peut-être on a donné le .bed pour basename
    if(length(grep("\\.bed$", basename)) > 0) {
      basename <- sub("\\.bed$", "", basename)
      bim <- paste(basename, ".bim", sep="")
      fam <- paste(basename, ".fam", sep="")
      bed <- path.expand(paste(basename, ".bed", sep=""))
      if(!is.null(rds)) rds <- paste(basename, ".rds", sep="")
    }
  }

  if(!file.exists(bed)) stop("file ", bed, " not found")

  if(is.null(rds) || !file.exists(rds)) {
    if(!file.exists(fam)) stop("file ", fam, " not found")
    if(!file.exists(bim)) stop("file ", bim, " not found")

    if(verbose) cat("Reading", fam, "\n")
    ped <- read.table(fam, stringsAsFactors = FALSE) 
    colnames(ped) <- pednames

    if(verbose) cat("Reading", bim, "\n")
    snp <- read.table(bim, stringsAsFactors = FALSE)
    colnames(snp) <- snpnames

    if(verbose) cat("Reading", bed, "\n")
    bed <- .Call(`_gaston_read_bed_file`, bed, nrow(ped), nrow(snp))
    x <- new("bed.matrix", bed = bed, snps = snp, ped = ped,                            
      p = NULL, mu = NULL, sigma = NULL, standardize_p = FALSE,
      standardize_mu_sigma = FALSE )
    if(getOption("gaston.auto.set.stats", TRUE)) x <- set.stats(x, verbose = verbose)
    return(x)
  }
 
  # reading rds [on ne fait pas set.stats dans ce cas !!]
  if(verbose) cat("Reading", rds, "\n")
  x <- readRDS(rds)
  if ( is(x) != "bed.matrix" ) stop("The object in file ", rds, " is not a bed.matrix")

  if(verbose) cat("Reading", bed, "\n")
  x@bed <- .Call(`_gaston_read_bed_file`, bed, nrow(x@ped), nrow(x@snps))

  x
}

