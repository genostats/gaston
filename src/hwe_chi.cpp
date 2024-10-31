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

