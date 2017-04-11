#include <Rcpp.h>
using namespace Rcpp;

// h2*s + (1-h2) > 0
// h2*(s-1) > -1
// si s > 1
// h2 > 1/(1-s) donc min_h2 = max 1/(1-s) pour s > 1
// si s < 1 h2 < 1/(1-s) donc max_h2 = min 1/(1-s) pour s < 1
// l'utilisateur fournit une valeur a priori et on met à jour en fonction
// de cette contrainte
// !!! valeur 1e-6 arbitraire pour que la vraisemblance reste bien définie aux bornes...
void min_max_h2(NumericVector Sigma, double & min_h2, double & max_h2) {
  int n = Sigma.size();
  // max_h2 = std::numeric_limits<double>::infinity();
  // min_h2 = -std::numeric_limits<double>::infinity();
  for(int i = 0; i < n; i++) {
    double s = Sigma[i];
    if(s > 1) {
      double u = 1/(1-s) + 1e-6;
      min_h2 = (min_h2 > u)?min_h2:u;
    }
    else if(s < 1) {
      double u = 1/(1-s) - 1e-6;
      max_h2 = (max_h2 < u)?max_h2:u;
    }
  }
}


