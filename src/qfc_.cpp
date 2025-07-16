#include <Rcpp.h>

extern "C" void qfc(double* lb1, double* nc1, int* n1, int *r1, double *sigma, double *c1, int *lim1, double *acc, double* trace, int* ifault, double *res);

// [[Rcpp::export]]
Rcpp::List qfc_(Rcpp::NumericVector lambda, Rcpp::NumericVector nonCentrality, Rcpp::IntegerVector degreeOfFreedom, double sigma, double c, int lim, double accuracy) {
  Rcpp::NumericVector trace(7);
  int ifault;
  double res;
  int r = lambda.size();
  if(r != nonCentrality.size() || r != degreeOfFreedom.size())
    Rcpp::stop("lambda, nonCentrality, and degreeOfFreedom should have same size\n");

  qfc(&lambda[0], &nonCentrality[0], &degreeOfFreedom[0], &r, &sigma, &c, &lim, &accuracy, &trace[0], &ifault, &res);

  Rcpp::List L;
  L["trace"] = trace;
  L["ifault"] = ifault;
  L["res"] = res;

  return L;
}


