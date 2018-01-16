#include <Rcpp.h>
using namespace Rcpp;

// lp = -log10(p) is assumed to be sorted in ascending order
//[[Rcpp::export]]
IntegerVector logp_thinning(NumericVector logp, double step) {
  std::vector<int> r;
  double x = 0;
  int i = 0, n = logp.size();
  debut:
    while(logp[i] < x && i < n) i++;
    if(i < n) {
      r.push_back( i+1 );
      x = logp[i] + step;
      goto debut;
    }
  return wrap(r);
}

RcppExport SEXP gg_logp_thinning(SEXP logpSEXP, SEXP stepSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type logp(logpSEXP);
    Rcpp::traits::input_parameter< double >::type step(stepSEXP);
    rcpp_result_gen = Rcpp::wrap(logp_thinning(logp, step));
    return rcpp_result_gen;
END_RCPP
}

