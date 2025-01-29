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

