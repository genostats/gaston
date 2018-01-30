
# doit marcher quelle que soit l'option choisie à la compilation
setThreadOptions <- function(numThreads = "auto", stackSize = "auto") {
  # 1 - on transmet à RcppParallel
  RcppParallel::setThreadOptions(numThreads, stackSize);
  # 2 - au cas où on a compilé sans RcppParallel (momentanément requis par le CRAN) on utilise notre mécanisme interne
  n <- numThreads
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
