#include <Rcpp.h>
using namespace Rcpp;

// lp = -log10(p) is assumed to be sorted in ascending order
//[[Rcpp::export]]
IntegerVector logp_thinning(NumericVector logp, double step) {
  std::vector<int> r;
  double x = 0;
  int i = 0, n = logp.size();
  debut:
    while(i < n && logp[i] < x)  i++;
    if(i < n) {
      if(std::isnan(logp[i])) 
        i++;
      else {
        x = logp[i] + step;
        r.push_back( ++i );
      }
      goto debut;
    }
  return wrap(r);
}

