.onLoad <- function(libname, pkgname) {
  options(gaston.chr.x = 23L, gaston.chr.y = 24L, gaston.chr.mt = 26L, gaston.autosomes = c(1:22,25L) )
  rnorm(1); # force seed initialisation (not done when the RNG is called from C++ code)
}
