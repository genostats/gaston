#include <Rcpp.h>
#include <RcppEigen.h>
#include "matrix-varia.h"

using namespace Rcpp;
typedef Map<MatrixXd> Map_MatrixXd;

double gg_pre_likelihood(NumericVector Y, NumericMatrix X, List K_, NumericVector h2) {
  Map_MatrixXd y(as<Map<MatrixXd> >(Y));
  Map_MatrixXd x(as<Map<MatrixXd> >(X));

  int s = K_.length();

  int n(y.rows()), p(x.cols());

  MatrixXd V(n,n), P(n,n), Vi(n,n);
  MatrixXd XViX(p,p), XViX_i(p,p);
  MatrixXd ViX(n,p);
  VectorXd Py(n);

  double log_detV, detV, d1, log_d1;

  // calcul de V0
  double residu = 1 - sum(h2);
  V = residu*MatrixXd::Identity(n,n);
  for(int j = 0; j < s; j++) V.noalias() += h2[j]*as<Map_MatrixXd>(as<NumericMatrix>(K_[j]));

  // Calcul de Vi = inverse(V0)
  sym_inverse(V,Vi,log_detV,detV,1e-7); // detruit V

  // Calcul de X' Vi X et de P0
  ViX.noalias() = Vi * x;
  XViX.noalias() = x.transpose() * ViX;
  sym_inverse(XViX, XViX_i, log_d1, d1, 1e-5);
  P.noalias() = Vi - ViX * XViX_i * ViX.transpose();

  // Py = P * y (en tenant ompte de la symmétrie de P)
  Py.noalias()   =  P.selfadjointView<Lower>() * y;
  double logL = -0.5*(log_detV + log_d1 + (n-p)*log(Py.dot(y.col(0))) + (n-p)*(1-log((double) (n-p))));
  return logL;
}

//[[Rcpp::export]]
double gg_pre_likelihood_nofix(NumericVector Y, List K_, NumericVector h2) {
  Map_MatrixXd y(as<Map<MatrixXd> >(Y));

  int s = K_.length();
  int n(y.rows());

  MatrixXd V(n,n), P(n,n), Vi(n,n);
  VectorXd Py(n);

  double log_detV, detV;

  // calcul de V 
  double residu = 1 - sum(h2);
  V = residu*MatrixXd::Identity(n,n);
  for(int j = 0; j < s; j++) V.noalias() += h2[j]*as<Map_MatrixXd>(as<NumericMatrix>(K_[j]));

  // Calcul de Vi = inverse(V)
  sym_inverse(V,Vi,log_detV,detV,1e-7); // detruit V

  // Py = P * y = Vi * y 
  Py.noalias()   =  Vi.selfadjointView<Lower>() * y;
  double logL = -0.5*(log_detV + n*log(Py.dot(y.col(0))) + n*(1-log((double)n)));
  return logL;

}

RcppExport SEXP gg_pre_likelihood(SEXP YSEXP, SEXP XSEXP, SEXP K_SEXP, SEXP h2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type Y(YSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< List >::type K_(K_SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type h2(h2SEXP);
    rcpp_result_gen = Rcpp::wrap(gg_pre_likelihood(Y, X, K_, h2));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP gg_pre_likelihood_nofix(SEXP YSEXP, SEXP K_SEXP, SEXP h2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type Y(YSEXP);
    Rcpp::traits::input_parameter< List >::type K_(K_SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type h2(h2SEXP);
    rcpp_result_gen = Rcpp::wrap(gg_pre_likelihood_nofix(Y, K_, h2));
    return rcpp_result_gen;
END_RCPP
}

