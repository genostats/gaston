#include <Rcpp.h>
#include <iostream>
#include <ctime>
#include <limits>

double hwe0(int a0, int a1, int a2) {
  int n = a0 + a1 + a2;
  if(n == 0)  // no data
    return std::numeric_limits<double>::quiet_NaN();
  if( (a0 == 0 && a1 == 0) || (a1 == 0 && a2 == 0)) // monomorphe
    return 1;
  if(a2 > a0) { // swap
    int tmp = a2;
    a2 = a0;
    a0 = tmp;
  }
  int m = 2*a2+a1; // allèles rares
  int s = ((2*n-m)*m)/(2*n);
  if(s%2 != m%2) s++;
  double grand_sum(1), small_sum(0);
  double target = 0; // juste pour supprimer le warning sur l'initiatilisation
  if(a1 < s) { // on calcule d'abord pour b1 < s
    int b2 = (m-s)/2;
    int b1 = s;
    int b0 = n - b1 - b2;
    double X = 1;
    bool over = false;
    while(b1 >= 2) {
      X *= (double) (b1*(b1-1)) / (4*(1+b0)*(1+b2));
      grand_sum += X;
      if(b1 == a1+2) {  
        target = X;
        over = true;
      }
      if(over) small_sum += X;
      b1 -= 2; b0++; b2++;
    }
    b2 = (m-s)/2;
    b1 = s;
    b0 = n - b1 - b2;
    X = 1;
    over = false;
    while(b1 <= m-2) {
      X *= (double) (4*(b0)*(b2)) / ((b1+2)*(b1+1)) ;
      grand_sum += X;
      if(over) 
        small_sum += X;
      else if(X <= target) {
        over = true;
        small_sum += X;
      }
      b1 += 2; b0--; b2--;
    }
  } else if(a1 > s) { // dans l'autre sens
    int b2 = (m-s)/2;
    int b1 = s;
    int b0 = n - b1 - b2;
    double X = 1;
    bool over = false; 
    while(b1 <= m-2) {
      X *= (double) (4*(b0)*(b2)) / ((b1+2)*(b1+1)) ;
      grand_sum += X;
      if(b1+2 == a1) {
        target = X;
        over = true;
      }
      if(over) small_sum += X;
      b1 += 2; b0--; b2--;
    }
    b2 = (m-s)/2;
    b1 = s;
    b0 = n - b1 - b2;
    X = 1;
    over = false;
    while(b1 >= 2) {
      X *= (double) (b1*(b1-1)) / (4*(1+b0)*(1+b2));
      grand_sum += X;
      if(over)
        small_sum += X;
      else if(X <= target) {
        over = true;
        small_sum += X;
      }
      b1 -= 2; b0++; b2++;
    }
  } else { // a1 = s!! -> target = 1
    int b2 = (m-s)/2;
    int b1 = s;
    int b0 = n - b1 - b2;
    double X = 1;
    bool over = false;
    target = 1;
    while(b1 <= m-2) {
      X *= (double) (4*(b0)*(b2)) / ((b1+2)*(b1+1)) ;
      grand_sum += X;
      if(over)
        small_sum += X;
      else if(X <= target) {
        over = true;
        small_sum += X;
      }
      b1 += 2; b0--; b2--;
    }
    b2 = (m-s)/2;
    b1 = s;
    b0 = n - b1 - b2;
    X = 1;
    over = false;
    while(b1 >= 2) {
      X *= (double) (b1*(b1-1)) / (4*(1+b0)*(1+b2));
      grand_sum += X;
      if(over)
        small_sum += X;
      else if(X <= target) {
        over = true;
        small_sum += X;
      }
      b1 -= 2; b0++; b2++;
    }
  }
  if(target >= 1) small_sum += 1; 
  return(small_sum/grand_sum);
}


// [[Rcpp::export]]
Rcpp::NumericVector hwe(Rcpp::IntegerVector N0, Rcpp::IntegerVector N1, Rcpp::IntegerVector N2) {
  int n = N0.size();
  Rcpp::NumericVector p(n);
  for(int i = 0; i < n; i++) p(i) = hwe0(N0(i), N1(i), N2(i));
  return p;
}

RcppExport SEXP gg_hwe(SEXP N0SEXP, SEXP N1SEXP, SEXP N2SEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type N0(N0SEXP );
        Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type N1(N1SEXP );
        Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type N2(N2SEXP );
        Rcpp::NumericVector __result = hwe(N0, N1, N2);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}


