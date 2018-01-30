.onAttach <- function(libname, pkgname) {
  options(gaston.chr.x = 23L, gaston.chr.y = 24L, gaston.chr.mt = 26L, gaston.autosomes = c(1:22,25L) )
  rnorm(1); # force RNG seed initialisation (not done when the RNG is called from C++ code)

  n <- RcppParallel::defaultNumThreads()
  n <- 2**max(0,floor(min(3, log2(n) - 1)))
  setThreadOptions(n)
  packageStartupMessage("Gaston set number of threads to ", n, ". Use setThreadOptions() to modify this.")
}

.onLoad <- function(libname, pkgname) {
  options(gaston.chr.x = 23L, gaston.chr.y = 24L, gaston.chr.mt = 26L, gaston.autosomes = c(1:22,25L) )
  rnorm(1); # force RNG seed initialisation (not done when the RNG is called from C++ code)
}

