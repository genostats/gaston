#include <Rcpp.h>
#include "matrix-varia.h"
#include "any.h"
#include "ai-reml-logit-nk.h"
#include "ai-reml-logit-nk-covar.h"

using namespace Rcpp;
typedef Map<MatrixXd> Map_MatrixXd;


//[[Rcpp::export]]
List AIREMLn_logit_nofix(NumericVector Y, List K_, bool constraint, NumericVector min_tau_, int max_iter, double eps, bool verbose,
                         NumericVector tau0, bool start_tau, bool get_P) {
						 
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
  
  AIREMLn_logit_nofix(y, K, constraint, min_tau, max_iter, eps, verbose, tau, niter, P, omega, start_tau);

  List L;

  L["tau"] = tau;
  L["niter"] = niter;
  if (get_P) L["P"] = P;
  L["BLUP_omega"] = omega;

  return L;
}

List AIREMLn_logit(NumericVector Y, NumericMatrix X, List K_, bool constraint, NumericVector min_tau_, int max_iter, double eps, bool verbose,
                   NumericVector tau0, NumericVector beta0, bool start_tau, bool start_beta, bool get_P) {
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

  AIREMLn_logit(y, x, K, constraint, min_tau, max_iter, eps, verbose, tau, niter, P, omega, beta, varbeta, start_tau, start_beta);

  List L;

  L["tau"] = tau;
  L["niter"] = niter;
  if (get_P) L["P"] = P;
  L["BLUP_omega"] = omega;
  L["BLUP_beta"] = beta;
  L["varbeta"] = varbeta;

  return L;
}


RcppExport SEXP gg_AIREMLn_logit_nofix(SEXP YSEXP, SEXP K_SEXP, SEXP constraintSEXP, SEXP min_tau_SEXP, SEXP max_iterSEXP, SEXP epsSEXP, SEXP verboseSEXP, SEXP tau0SEXP, SEXP start_tauSEXP, SEXP get_PSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericVector >::type Y(YSEXP );
        Rcpp::traits::input_parameter< List >::type K_(K_SEXP );
        Rcpp::traits::input_parameter< bool >::type constraint(constraintSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type min_tau_(min_tau_SEXP );
        Rcpp::traits::input_parameter< int >::type max_iter(max_iterSEXP );
        Rcpp::traits::input_parameter< double >::type eps(epsSEXP );
        Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP );
		Rcpp::traits::input_parameter< NumericVector >::type tau0(tau0SEXP );
        Rcpp::traits::input_parameter< bool >::type start_tau(start_tauSEXP );
	Rcpp::traits::input_parameter< bool >::type get_P( get_PSEXP );
        List __result = AIREMLn_logit_nofix(Y, K_, constraint, min_tau_, max_iter, eps, verbose, tau0, start_tau, get_P);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}

RcppExport SEXP gg_AIREMLn_logit(SEXP YSEXP, SEXP XSEXP, SEXP K_SEXP, SEXP constraintSEXP, SEXP min_tau_SEXP, SEXP max_iterSEXP, SEXP epsSEXP, SEXP verboseSEXP, SEXP tau0SEXP, SEXP beta0SEXP, SEXP start_tauSEXP, SEXP start_betaSEXP, SEXP get_PSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericVector >::type Y(YSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP );
        Rcpp::traits::input_parameter< List >::type K_(K_SEXP );
        Rcpp::traits::input_parameter< bool >::type constraint(constraintSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type min_tau_(min_tau_SEXP );
        Rcpp::traits::input_parameter< int >::type max_iter(max_iterSEXP );
        Rcpp::traits::input_parameter< double >::type eps(epsSEXP );
        Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP );
		Rcpp::traits::input_parameter< NumericVector >::type tau0(tau0SEXP );
        Rcpp::traits::input_parameter< NumericVector >::type beta0(beta0SEXP );
        Rcpp::traits::input_parameter< bool >::type start_tau(start_tauSEXP );
        Rcpp::traits::input_parameter< bool >::type start_beta(start_betaSEXP );
	Rcpp::traits::input_parameter< bool >::type get_P(get_PSEXP );
        List __result = AIREMLn_logit(Y, X, K_, constraint, min_tau_, max_iter, eps, verbose, tau0, beta0, start_tau, start_beta, get_P);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}

