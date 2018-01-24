
# doit marcher quelle que soit l'option choisie à la compilation
setThreadOptions <- function(...) {
  # 1 - on transmet à RcppParallel
  RcppParallel::setThreadOptions(...);
  # 2 - au cas où on a compilé sans RcppParallel (momentanément requis par le CRAN) on utilise notre mécanisme interne
  n <- list(...)
  if(is.null(names(n))) 
    n <- n[[1]]
  else if("numThreads" %in% names(n))
    n <- n$numThreads
  else if(names(n)[1] == "")
    n <- n[[1]]

  if(n == "auto") {
    .Call("set_nb_threads_to_default", PACKAGE = "gaston")
  } else {
    n <- as.integer(n)[1]
    if(is.na(n)) 
      warning("Bad argument")
    else
      .Call("set_nb_threads",  PACKAGE = "gaston", n)
  }
}
