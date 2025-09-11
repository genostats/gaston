#' Set Genetic Distance
#' 
#' @description
#' Returns an updated \code{\link{bed.matrix}} with genetic
#' distances in centimorgan computed from the variant positions
#' 
#' @param x  A \code{\link{bed.matrix}}
#' @param map  The genetic map, given by a list of data frames (see Details)
#' @param verbose  If \code{TRUE}, display information on the function actions
#' 
#' 
#' @details
#' A map is a list of data frames, with names corresponding to chromosomes.
#' Each of these data frames must have columns \code{pos} and \code{dist} corresponding
#' to positions in bp and cM, respectively.
#' 
#' Such maps are too large to be included in a CRAN package. You can get two genetic
#' maps for the Human Genome (build 36 and 37) in the package \code{HumanGeneticMap}
#' on GitHub.
#' 
#' To install this package, run
#' 
#' \code{ install.packages("HumanGeneticMap", repos="https://genostats.github.io/R/") }
#' 
#' You can then use this function with \code{set.dist(x, HumanGeneticMap::genetic.map.b36)}
#' for example, for positions on the build 36. Use \code{map = HumanGeneticMap::genetic.map.b37})
#' for the build 37 and \code{map = HumanGeneticMap::genetic.map.b38}) for the build 38.
#' 
#' 
#' @return
#' A \code{\link{bed.matrix}} similar to \code{x}, with updated values in \code{x@snps$dist}.
#' 
#' 
#' @export set.dist
set.dist <- function(x, map, verbose =  getOption("gaston.verbose",TRUE)) {
  if(missing(map) || !is.list(map)) 
    stop("This functions needs a genetic map. Please read the manual.")
  for(chr in as.character(unique(x@snps$chr))) {
    MAP <- map[[chr]]
    if( is.null(MAP) ) {
      if(verbose) cat("No map for chromosome", chr, "\n")
      next
    }
    if(verbose) cat("Setting distances for chromosome", chr, "\n")
    w <- which(x@snps$chr == chr)
    f <- approxfun( MAP$pos, MAP$dist, rule = 2)
    x@snps$dist[w] <- f( x@snps$pos[w] )
  }
  x
}
