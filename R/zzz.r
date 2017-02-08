.onLoad <- function(libname, pkgname) {
  options(gaston.chr.x = 23, gaston.chr.y = 24, gaston.chr.mt = 26, gaston.autosomes = c(1:22,25) )
  rnorm(1); # force seed initialisation (not done when the RNG is called from C++ code)
}
