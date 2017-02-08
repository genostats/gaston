#include <Rcpp.h>
#include <iostream>
#include "matrix4.h"

using namespace Rcpp;

// [[Rcpp::export]]
XPtr<matrix4> bind_snps(List L) {
  int s = L.size();
  if(s < 2) Rf_error("Can't bind less than two matrices!");
  XPtr<matrix4> first = as<XPtr<matrix4> >(L[0]);
  int n = first->ncol;
  int m = first->nrow;
  for(int i = 1; i < s; i++) {
    XPtr<matrix4> nxt = as<XPtr<matrix4> >(L[i]);
    if(n != nxt->ncol) Rf_error("Dimensions mismatch");
    m += nxt->nrow;
  }
  XPtr<matrix4> r(new matrix4(m,n));
  int k = 0;
  for(int i = 0; i < s; i++) {
    XPtr<matrix4> nxt = as<XPtr<matrix4> >(L[i]);
    for(int ii = 0; ii < nxt->nrow; ii++) {
      std::copy(nxt->data[ii], nxt->data[ii]+nxt->true_ncol, r->data[k++]);
    }
  }
  return r;
}

RcppExport SEXP gg_bind_snps(SEXP LSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< List >::type L(LSEXP );
        XPtr<matrix4> __result = bind_snps(L);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}

