#include <Rcpp.h>
#include <RcppEigen.h>
#include "matrix-varia.h"

using namespace Rcpp;
typedef Map<MatrixXd> Map_MatrixXd;

double gg_re_likelihood(NumericVector Y, NumericMatrix X, List K_, NumericVector theta) {
  Map_MatrixXd y(as<Map<MatrixXd> >(Y));
  Map_MatrixXd x(as<Map<MatrixXd> >(X));

  int s = K_.length();

  /*  std::vector<Map_MatrixXd> K;
  for(int i = 0; i < s; i++)
    K.push_back(as<Map_MatrixXd>(as<NumericMatrix>(K_[i])));

  return re_likelihood(y, x, K, theta);
  */
  int n(y.rows()), p(x.cols());

  MatrixXd V(n,n), P(n,n), Vi(n,n);
  MatrixXd XViX(p,p), XViX_i(p,p);
  MatrixXd ViX(n,p);
  VectorXd Py(n);

  double log_detV, detV, d1, log_d1;

  // calcul de V 
  V = theta(0)*MatrixXd::Identity(n,n);
  for(int j = 0; j < s; j++) V.noalias() += theta[j+1]*as<Map_MatrixXd>(as<NumericMatrix>(K_[j]));

  // Calcul de Vi = inverse(V)
  sym_inverse(V,Vi,log_detV,detV,1e-7); // detruit V

  // Calcul de X' Vi X et de P
  ViX.noalias() = Vi * x;
  XViX.noalias() = x.transpose() * ViX;
  sym_inverse(XViX, XViX_i, log_d1, d1, 1e-5);
  P.noalias() = Vi - ViX * XViX_i * ViX.transpose();

  // Py = P * y (en tenant ompte de la symmétrie de P)
  Py.noalias()   =  P.selfadjointView<Lower>() * y;
  double logL = -0.5*(log_detV + log_d1 + Py.dot(y.col(0)));
  return logL;
}

//[[Rcpp::export]]
double gg_re_likelihood_nofix(NumericVector Y, List K_, NumericVector theta) {
  Map_MatrixXd y(as<Map<MatrixXd> >(Y));

  int s = K_.length();
  /*  std::vector<Map_MatrixXd> K;
  for(int i = 0; i < s; i++)
    K.push_back(as<Map_MatrixXd>(as<NumericMatrix>(K_[i])));

  return re_likelihood_nofix(y, K, theta);
  */

  int n(y.rows());

  MatrixXd V(n,n), P(n,n), Vi(n,n);
  VectorXd Py(n);

  double log_detV, detV;

  // calcul de V 
  V = theta(0)*MatrixXd::Identity(n,n);
  for(int j = 0; j < s; j++) V.noalias() += theta[j+1]*as<Map_MatrixXd>(as<NumericMatrix>(K_[j]));

  // Calcul de Vi = inverse(V)
  sym_inverse(V,Vi,log_detV,detV,1e-7); // detruit V

  // Py = P * y = Vi * y 
  Py.noalias()   =  Vi.selfadjointView<Lower>() * y;
  double logL = -0.5*(log_detV + Py.dot(y.col(0)));
  return logL;

}

RcppExport SEXP gg_re_likelihood(SEXP YSEXP, SEXP XSEXP, SEXP K_SEXP, SEXP thetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type Y(YSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< List >::type K_(K_SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type theta(thetaSEXP);
    rcpp_result_gen = Rcpp::wrap(gg_re_likelihood(Y, X, K_, theta));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP gg_re_likelihood_nofix(SEXP YSEXP, SEXP K_SEXP, SEXP thetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type Y(YSEXP);
    Rcpp::traits::input_parameter< List >::type K_(K_SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type theta(thetaSEXP);
    rcpp_result_gen = Rcpp::wrap(gg_re_likelihood_nofix(Y, K_, theta));
    return rcpp_result_gen;
END_RCPP
}

