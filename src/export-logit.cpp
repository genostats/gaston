#include <Rcpp.h>
#include "matrix-varia.h"
#include "logit_model.h"

using namespace Rcpp;
typedef Map<MatrixXd> Map_MatrixXd;

// [[Rcpp::export]]
List logistic(NumericVector Y, NumericMatrix X, double eps) {
  Map_MatrixXd y(as<Map<MatrixXd> >(Y));
  Map_MatrixXd x(as<Map<MatrixXd> >(X));
  
  VectorXd beta(x.cols());
  MatrixXd varbeta(x.cols(),x.cols());

  beta.setZero();
  logistic_model2<double>(y, x, beta, varbeta, eps);

  List L;
  L["beta"] = beta;
  L["varbeta"] = varbeta;

  return L;
}

RcppExport SEXP gg_logistic(SEXP YSEXP, SEXP XSEXP, SEXP epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type Y(YSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    rcpp_result_gen = Rcpp::wrap(logistic(Y, X, eps));
    return rcpp_result_gen;
END_RCPP
}
