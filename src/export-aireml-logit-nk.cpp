#include <Rcpp.h>
#include "matrix-varia.h"
#include "any.h"
#include "ai-reml-logit-nk.h"
#include "ai-reml-logit-nk-covar.h"

using namespace Rcpp;
typedef Map<MatrixXd> Map_MatrixXd;


//[[Rcpp::export]]
List AIREMLn_logit_nofix(NumericVector Y, List K_, bool constraint, NumericVector min_tau_, int max_iter, double eps, bool verbose,
                         NumericVector tau0, bool start_tau, bool get_P, bool EM) {

  Map_MatrixXd y(as<Map_MatrixXd>(Y));

  int s = K_.length();
  std::vector<Map_MatrixXd> K;
  for(int i = 0; i < s; i++) 
    K.push_back(as<Map_MatrixXd>(as<NumericMatrix>(K_[i])));

  int niter;
  VectorXd tau(s), omega(y.rows());
  MatrixXd P(y.rows(),y.rows());
  Map_MatrixXd min_tau(as<Map_MatrixXd>(min_tau_));
  
  for(int i=0; i < s; i++) tau(i) = tau0(i);
  
  AIREMLn_logit_nofix(y, K, constraint, min_tau, max_iter, eps, verbose, tau, niter, P, omega, start_tau, EM);

  List L;

  L["tau"] = tau;
  L["niter"] = niter;
  if (get_P) L["P"] = P;
  L["BLUP_omega"] = omega;

  return L;
}

//[[Rcpp::export]]
List AIREMLn_logit(NumericVector Y, NumericMatrix X, List K_, bool constraint, NumericVector min_tau_, int max_iter, double eps, bool verbose,
                   NumericVector tau0, NumericVector beta0, bool start_tau, bool start_beta, bool get_P, bool EM) {

  Map_MatrixXd y(as<Map<MatrixXd> >(Y));
  Map_MatrixXd x(as<Map<MatrixXd> >(X));

  int s = K_.length();
  std::vector<Map_MatrixXd> K;
  for(int i = 0; i < s; i++)
    K.push_back(as<Map_MatrixXd>(as<NumericMatrix>(K_[i])));

  int niter;
  VectorXd tau(s), omega(y.rows()), beta(x.cols());
  MatrixXd varbeta(x.cols(),x.cols()), P(y.rows(),y.rows());
  Map_MatrixXd min_tau(as<Map_MatrixXd>(min_tau_));
  
  for(int i=0; i < s; i++) tau(i) = tau0(i);
  for(int i=0; i < x.cols(); i++) beta(i) = beta0(i);

  AIREMLn_logit(y, x, K, constraint, min_tau, max_iter, eps, verbose, tau, niter, P, omega, beta, varbeta, start_tau, start_beta, EM);

  List L;

  L["tau"] = tau;
  L["niter"] = niter;
  if (get_P) L["P"] = P;
  L["BLUP_omega"] = omega;
  L["BLUP_beta"] = beta;
  L["varbeta"] = varbeta;

  return L;
}

