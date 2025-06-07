#' Save a \code{\link{bed.matrix}}
#' 
#' @description
#' Save a \code{\link{bed.matrix}} in several files
#' 
#' @param x  A \code{\link{bed.matrix}}
#' @param basename  Basename of all files
#' @param bed  Name of the \code{.bed} file
#' @param fam  Name of the \code{.fam} file
#' @param bim  Name of the \code{.bim} file
#' @param rds  Name of the \code{.rds} file
#' 
#' 
#' @details  If any of \code{bed}, \code{fam}, \code{bim} and \code{rds} is \code{NULL},
#' the corresponding file will not be written.
#' 
#' The \code{.fam} and \code{.bim} files are useful for reading files with other softwares.
#' The \code{.rds} file can be read by \code{read.bed.matrix}.
#' 
#' The \code{.bed}, \code{.fam} and \code{.bim} files follow the PLINK specifications
#' (\url{http://zzz.bwh.harvard.edu/plink/binary.shtml}).
#' 
#' 
#' 
#' @seealso  \code{\link{read.bed.matrix}}, \code{\link[base:saveRDS]{saveRDS}}
#' 
#' 
#' @examples
#' 
#' # Load data
#' data(LCT)
#' x <- as.bed.matrix(LCT.gen, LCT.fam, LCT.bim)
#' 
#' # Write object in LCT.bed and LCT.RData
#' \dontrun{
#' write.bed.matrix(x, "LCT")
#' }
#' 
#' @export write.bed.matrix
write.bed.matrix <- function(x, basename, bed = paste(basename, ".bed", sep=""), fam = paste(basename, ".fam", sep=""),
                                          bim = paste(basename, ".bim", sep=""), rds = paste(basename, ".rds", sep="")) {
  if ( is(x) != "bed.matrix" ) stop("x must be a bed.matrix")

  if(!is.null(fam)) {
    if(all(pednames %in% names(x@ped)))
      write.table(x@ped[,pednames], file = fam, row.names=FALSE, col.names=FALSE, quote = FALSE)
    else
      warning("Can't create .fam file from x: missing columns in x@ped")
  }

  if(!is.null(bim)){
    if(all(snpnames %in% names(x@snps)))
      write.table(x@snps[,snpnames], file = bim, row.names=FALSE, col.names=FALSE, sep = "\t", quote = FALSE)
    else
      warning("Can't create .bim file from x: missing columns in x@snps")
  }

  if(!is.null(rds))
    saveRDS(x, rds)

  if(!is.null(bed)) {
    bed <- path.expand(bed)
    invisible(.Call(`_gaston_write_bed_file`, x@bed, bed))
  }
}

