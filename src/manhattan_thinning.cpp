#include <Rcpp.h>
#include <limits>

using namespace Rcpp;

double min_(NumericVector x) { // suite à problèmes avec na_omit qui ne reconnaît pas les -NA_real_ ... (Rcpp 0.12.14)
  double m = std::numeric_limits<double>::infinity();
  int n = x.length();
  for(int i = 0; i < n; i++) {
    if(x[i] < m) // false si x[i] est nan ..
      m = x[i];
  }
  return m;
}

double max_(NumericVector x) {
  double m = -std::numeric_limits<double>::infinity();
  int n = x.length();
  for(int i = 1; i < n; i++) {
    if(x[i] > m) 
      m = x[i];
  }
  return m;
}

// les valeurs sont dans l'ordre où elles vont être plottés
// ie coord est trié dans l'ordre croissant
// [[Rcpp::export]]
IntegerVector manhattan_thinning(NumericVector x, NumericVector y, int mx, int my) {
  int n = x.length();
  if(n != y.length()) stop("x and y should have the same length");
  double min_y = min_(y);
  double dx = (max_(x) - min_(x))/mx;
  double dy = (max_(y) - min_y)/my;
  std::vector<int> Y(my+1);
  
  std::vector<int> r;
  for(int i = 0; i < n; i++) {
    if(std::isnan(x[i]) || std::isnan(y[i])) 
      continue;
    int j = (y[i] - min_y)/dy;
    if( Y[j] == 0 || (x[i] - x[Y[j]-1]) > dx ) {
      r.push_back(i+1);
      Y[j] = i+1;
    }
  }
  return wrap(r);
}

RcppExport SEXP gg_manhattan_thinning(SEXP xSEXP, SEXP ySEXP, SEXP mxSEXP, SEXP mySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< int >::type mx(mxSEXP);
    Rcpp::traits::input_parameter< int >::type my(mySEXP);
    rcpp_result_gen = Rcpp::wrap(manhattan_thinning(x, y, mx, my));
    return rcpp_result_gen;
END_RCPP
}

