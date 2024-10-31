#include <Rcpp.h>
#include <limits>

using namespace Rcpp;

double min_(NumericVector x) { // suite à problèmes avec na_omit qui ne reconnaît pas les -NA_real_ ... (Rcpp 0.12.14)
  double infinity = std::numeric_limits<double>::infinity();
  double m = infinity;
  int n = x.length();
  for(int i = 0; i < n; i++) {
    if(x[i] < m && x[i] > -infinity) // false si x[i] est nan .. [ +éviter les range infinis]
      m = x[i];
  }
  return m;
}

double max_(NumericVector x) {
  double infinity = std::numeric_limits<double>::infinity();
  double m = -infinity;
  int n = x.length();
  for(int i = 1; i < n; i++) {
    if(x[i] > m && x[i] < infinity) // idem
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
    if( std::isnan(x[i]) || std::isinf(x[i]) || std::isnan(y[i]) || std::isinf(y[i]) ) 
      continue;
    int j = (y[i] - min_y)/dy;
    if( Y[j] == 0 || (x[i] - x[Y[j]-1]) > dx ) {
      r.push_back(i+1);
      Y[j] = i+1;
    }
  }
  return wrap(r);
}
