#include <Rcpp.h>
#include <iostream>
#include "matrix4.h"

using namespace Rcpp;

inline uint8_t flip_strand2(uint8_t x) {
  if(x==3) return 3;
  return(2-x);
}

//[[Rcpp::export]]
XPtr<matrix4> bind_inds2(List L, LogicalMatrix flip) {
  int s = L.size();
  if(s < 2) 
    stop("Can't bind less than two matrices!");
  if(flip.nrow() != s)
    stop("Dimensions mismatch");

  XPtr<matrix4> first = as<XPtr<matrix4> >(L[0]);
  int n = first->ncol;
  int m = first->nrow;
  if(flip.ncol() != m)
    stop("Dimensions mismatch");
  for(int i = 1; i < s; i++) {
    XPtr<matrix4> nxt = as<XPtr<matrix4> >(L[i]);
    if(m != nxt->nrow) stop("Dimensions mismatch");
    n += nxt->ncol;
  }
  XPtr<matrix4> r(new matrix4(m,n));
  for(int i = 0; i < m; i++) {
    int k = 0;
    for(int j = 0; j < s; j++) {
      XPtr<matrix4> nxt = as<XPtr<matrix4> >(L[j]);
      for(int jj = 0; jj < nxt->ncol; jj++) {
        if(LogicalVector::is_na(flip(j,i))) {
          (*r)(i,k++) = 3; // NA
        } else if(flip(j,i)) {
          (*r)(i,k++) = flip_strand2( (*nxt)(i,jj) );
        } else {
          (*r)(i,k++) = (*nxt)(i,jj);
        }
      }
    }
  }
  return r;
}



RcppExport SEXP gg_bind_inds2(SEXP LSEXP, SEXP flipSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< List >::type L(LSEXP);
    Rcpp::traits::input_parameter< LogicalMatrix >::type flip(flipSEXP);
    __result = Rcpp::wrap(bind_inds2(L, flip));
    return __result;
END_RCPP
}

