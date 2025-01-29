logp.thinning <- function( lp, step ) {
   .Call(`_gaston_logp_thinning`, PACKAGE = "gaston", lp, step);
}
