
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


