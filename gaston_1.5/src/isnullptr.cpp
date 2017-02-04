#include <Rcpp.h>
#include <Rinternals.h>

RcppExport SEXP isnullptr(SEXP pointer) {
  SEXP result;
  bool x = (bool) R_ExternalPtrAddr(pointer);
  PROTECT( result = Rcpp::wrap(!x) );
  UNPROTECT(1); // tout ça est-il bien nécessaire ?!
  return result;
}

