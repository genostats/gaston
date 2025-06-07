#' Genomic Sex
#' 
#' @description
#' Returns an updated \code{\link{bed.matrix}} with a new variable for the genomic sex
#' of each individual.
#' 
#' @param x  A \code{\link{bed.matrix}}
#' @param plot  If \code{TRUE}, plots the variables used for the clustering
#' @param verbose  If \code{TRUE}, displays information on the function actions
#' 
#' @details
#' For each individual, the function uses the hetorozygosity rate for SNPs on X chromosome,
#' and the call rate for SNPs on the Y chromosomes (both statistics computed by \code{\link{set.stats}}),
#' to cluster the individuals using \code{\link{kmeans}}.
#' 
#' If \code{plot = TRUE}, a plot is produced with the two variables used and the clusters
#' determined by \code{\link{kmeans}}.
#' 
#' @return
#' A \code{\link{bed.matrix}} similar to \code{x}, with a new variable \code{x@ped$genomic.sex}
#' containing the genomic sex for each individual.
#' @seealso  \code{\link{set.stats}}, \code{\link{set.hwe}}
#' 
#' @keywords  Genomic Sex
#' @export set.genomic.sex
set.genomic.sex <- function(x, plot = FALSE, verbose = getOption("gaston.verbose",TRUE)) { 
  if( !all(c("hz.x", "callrate.y") %in% names(x@ped) )) {
    if(verbose) cat("Computing individual X heterozygosity and Y callrate stats\n")
    x <- set.stats.ped(x)
  }

  y.callrate <- x@ped$callrate.y
  x.hz       <- x@ped$hz.x

  # centres naturels des clusters
  # H : x.hz = 0 y.callrate = 1
  # F : x.hz = ? y.callrate = 0 
  max.x.hz <- max(x.hz)

  cluster <- kmeans( x = cbind( c(y.callrate,1,0) , c(x.hz,0,max.x.hz) ), 
                     centers = matrix(c(1, 0, 0, max.x.hz), ncol = 2) )$cluster
  x@ped$genomic.sex <- cluster[ 1:nrow(x) ]
  if(plot) {
    plot( x.hz, y.callrate, xlab = "X heterozygosity", ylab = "Y callrate", 
          col = x@ped$genomic.sex, main = "Determination of Genomic Sex", xlim=c(0,1), ylim=c(0,1), asp=1 )
    legend("bottomleft", pch = 1, col = 1:2, c("M","F"))
  }
  if(verbose) {
    cat("Slot @ped$genomic.sex has been set\n")
    diff.mf <- sum(x@ped$sex == 1 & x@ped$genomic.sex == 2)
    diff.fm <- sum(x@ped$sex == 2 & x@ped$genomic.sex == 1)
    sex.undef <- sum(x@ped$sex != 1 & x@ped$sex != 2)
    if(diff.mf > 0)
      cat("Found", diff.mf, "individuals with sex = M and genomic sex = F\n")
    if(diff.fm > 0)
      cat("Found", diff.fm, "individuals with sex = F and genomic sex = M\n")
    if(sex.undef > 0)
      cat("Found", sex.undef, "individuals with undefined sex\n");
  }
  x
}


