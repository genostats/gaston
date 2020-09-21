#include <Rcpp.h>
#include "matrix-varia.h"
#include "ai-reml-logit-1k.h"
#include "ai-reml-logit-1k-covar.h"

using namespace Rcpp;
typedef Map<MatrixXd> Map_MatrixXd;

// [[Rcpp::export]]
List AIREML1_logit_nofix(NumericVector Y, NumericMatrix K_, bool constraint, double min_tau, int max_iter, double eps, bool verbose,
                         double tau0, bool start_tau, bool get_P, bool EM) {
  Map_MatrixXd y(as<Map<MatrixXd> >(Y));
  Map_MatrixXd K(as<Map<MatrixXd> >(K_));

  int niter;
  MatrixXd P(y.rows(),y.rows());
  VectorXd omega(y.rows());
  double tau(tau0);

  AIREML1_logit_nofix(y, K, constraint, min_tau, max_iter, eps, verbose, tau, niter, P, omega, start_tau, EM);

  List L;
  L["tau"] = tau;
  L["niter"] = niter;
  if(get_P) L["P"] = P; // réalise la copie 
  L["BLUP_omega"] = omega;

  return L;
}

// [[Rcpp::export]]
List AIREML1_logit(NumericVector Y, NumericMatrix X, NumericMatrix K_, bool constraint, double min_tau, int max_iter, double eps, bool verbose,
                   double tau0, NumericVector beta0, bool start_tau, bool start_beta, bool get_P, bool EM) {
  Map_MatrixXd y(as<Map<MatrixXd> >(Y));
  Map_MatrixXd x(as<Map<MatrixXd> >(X));
  Map_MatrixXd K(as<Map<MatrixXd> >(K_));
  
  int niter;
  MatrixXd P(y.rows(),y.rows());
  VectorXd omega(y.rows());
  VectorXd beta(x.cols());
  for (int j=0; j<x.cols(); j++) beta(j)=beta0(j);
  double tau(tau0);
  MatrixXd varbeta(x.cols(),x.cols());

  AIREML1_logit(y, x, K, constraint, min_tau, max_iter, eps, verbose, tau, niter, P, omega, beta, varbeta, start_tau, start_beta, EM);

  List L;
  L["tau"] = tau;
  L["niter"] = niter;
  if(get_P) L["P"] = P;
  L["BLUP_omega"] = omega;
  L["BLUP_beta"] = beta;
  L["varbeta"] = varbeta;

  return L;
}




// AIREML1_logit_nofix
RcppExport SEXP gg_AIREML1_logit_nofix(SEXP YSEXP, SEXP K_SEXP, SEXP constraintSEXP, SEXP min_tauSEXP, SEXP max_iterSEXP, SEXP epsSEXP, SEXP verboseSEXP, SEXP tau0SEXP, SEXP start_tauSEXP, SEXP get_PSEXP, SEXP EMSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type Y(YSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type K_(K_SEXP);
    Rcpp::traits::input_parameter< bool >::type constraint(constraintSEXP);
    Rcpp::traits::input_parameter< double >::type min_tau(min_tauSEXP);
    Rcpp::traits::input_parameter< int >::type max_iter(max_iterSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< double >::type tau0(tau0SEXP);
    Rcpp::traits::input_parameter< bool >::type start_tau(start_tauSEXP);
    Rcpp::traits::input_parameter< bool >::type get_P(get_PSEXP);
    Rcpp::traits::input_parameter< bool >::type EM(EMSEXP);
    rcpp_result_gen = Rcpp::wrap(AIREML1_logit_nofix(Y, K_, constraint, min_tau, max_iter, eps, verbose, tau0, start_tau, get_P, EM));
    return rcpp_result_gen;
END_RCPP
}

// AIREML1_logit
RcppExport SEXP gg_AIREML1_logit(SEXP YSEXP, SEXP XSEXP, SEXP K_SEXP, SEXP constraintSEXP, SEXP min_tauSEXP, SEXP max_iterSEXP, SEXP epsSEXP, SEXP verboseSEXP, SEXP tau0SEXP, SEXP beta0SEXP, SEXP start_tauSEXP, SEXP start_betaSEXP, SEXP get_PSEXP, SEXP EMSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type Y(YSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type K_(K_SEXP);
    Rcpp::traits::input_parameter< bool >::type constraint(constraintSEXP);
    Rcpp::traits::input_parameter< double >::type min_tau(min_tauSEXP);
    Rcpp::traits::input_parameter< int >::type max_iter(max_iterSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< double >::type tau0(tau0SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type beta0(beta0SEXP);
    Rcpp::traits::input_parameter< bool >::type start_tau(start_tauSEXP);
    Rcpp::traits::input_parameter< bool >::type start_beta(start_betaSEXP);
    Rcpp::traits::input_parameter< bool >::type get_P(get_PSEXP);
    Rcpp::traits::input_parameter< bool >::type EM(EMSEXP);
    rcpp_result_gen = Rcpp::wrap(AIREML1_logit(Y, X, K_, constraint, min_tau, max_iter, eps, verbose, tau0, beta0, start_tau, start_beta, get_P, EM));
    return rcpp_result_gen;
END_RCPP
}

