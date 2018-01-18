manhattan.thinning <- function(x, y, mx, my ) {
  if( is.unsorted(x) ) {
    o <- order(x);
    x1 <- x[o]
    y1 <- y[o]
    w <- .Call('gg_manhattan_thinning', PACKAGE = "gaston", x1, y1, mx, my)
    return(o[w]);
  } 
  return(.Call('gg_manhattan_thinning', PACKAGE = "gaston", x, y, mx, my))
}
