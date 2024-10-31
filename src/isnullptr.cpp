#include <Rcpp.h>
#include <Rinternals.h>

// [[Rcpp::export]]
bool isnullptr(SEXP pointer) {
  bool x = (bool) R_ExternalPtrAddr(pointer);
  return !x;
}

