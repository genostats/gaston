# x = une bed.matrix
# map 
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
