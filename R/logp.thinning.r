logp.thinning <- function( lp, step ) {
   .Call('gg_logp_thinning', PACKAGE = "gaston", lp, step);
}
