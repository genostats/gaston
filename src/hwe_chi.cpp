#include <Rcpp.h>
#include <iostream>
#include <ctime>
#include <limits>

double hwe_chi0(unsigned int a0, unsigned int a1, unsigned int a2) {
  int n = a0 + a1 + a2;
  if(n == 0) // no data
    return std::numeric_limits<double>::quiet_NaN(); 
  double p = (double) (2*a2+a1) / (2*n);
  if(p == 0 || p == 1) return 1; // monomorphe
  double e0 = n*(1-p)*(1-p), e1 = 2*n*p*(1-p), e2 = n*p*p;
  double chi = (a0-e0)*(a0-e0)/e0 + (a1-e1)*(a1-e1)/e1 + (a2-e2)*(a2-e2)/e2;
  return R::pchisq(chi, 1., 0, 0);
}


// [[Rcpp::export]]
Rcpp::NumericVector hwe_chi(Rcpp::IntegerVector N0, Rcpp::IntegerVector N1, Rcpp::IntegerVector N2) {
  int n = N0.size();
  Rcpp::NumericVector p(n);
  for(int i = 0; i < n; i++) p(i) = hwe_chi0(N0(i), N1(i), N2(i));
  return p;
}

RcppExport SEXP gg_hwe_chi(SEXP N0SEXP, SEXP N1SEXP, SEXP N2SEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type N0(N0SEXP );
        Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type N1(N1SEXP );
        Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type N2(N2SEXP );
        Rcpp::NumericVector __result = hwe_chi(N0, N1, N2);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}


