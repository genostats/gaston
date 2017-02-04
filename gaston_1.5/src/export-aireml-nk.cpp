#include <Rcpp.h>
#include "matrix-varia.h"
#include "any.h"
#include "ai-reml-nk.h"
#include "ai-reml-nk-covar.h"

using namespace Rcpp;
typedef Map<MatrixXd> Map_MatrixXd;

//[[Rcpp::export]]
List AIREMLn_nofix(NumericVector Y, List K_, int EMsteps, int EMsteps_fail, double EM_alpha, bool constraint, 
                   double min_s2, NumericVector min_tau_, int max_iter, double eps, bool verbose,
                   NumericVector theta_, bool start_theta, bool get_P) {
  Map_MatrixXd y(as<Map_MatrixXd>(Y));

  int s = K_.length();
  std::vector<Map_MatrixXd> K;
  for(int i = 0; i < s; i++) 
    K.push_back(as<Map_MatrixXd>(as<NumericMatrix>(K_[i])));

  double logL, logL0, normgrad;
  int niter;
  MatrixXd P(y.rows(),y.rows());
  VectorXd theta(s+1), Py(y.rows()), omega(y.rows());
  Map_MatrixXd min_tau(as<Map_MatrixXd>(min_tau_));

  for(int i=0; i < s+1; i++) theta(i) = theta_(i);

  AIREMLn_nofix(y, K, EMsteps, EMsteps_fail, EM_alpha, constraint, min_s2, min_tau, max_iter, eps, 
                verbose, theta, logL, logL0, niter, normgrad, P, Py, omega, start_theta);

  List L;

  L["sigma2"] = theta(0);
  L["tau"] = theta.tail(s);
  L["logL"] = logL;
  L["logL0"] = logL0;
  L["niter"] = niter;
  L["norm_grad"] = normgrad;
  if(get_P) L["P"] = P;  
  L["Py"] = Py;  
  L["BLUP_omega"] = omega;

  return L;
}

//[[Rcpp::export]]
List AIREMLn_contrast(NumericVector Y, NumericMatrix X, List K_, int EMsteps, int EMsteps_fail, double EM_alpha, bool constraint,
             double min_s2, NumericVector min_tau_, int max_iter, double eps, bool verbose,
             NumericVector theta_, bool start_theta, bool get_P) {
  Map_MatrixXd y(as<Map<MatrixXd> >(Y));
  Map_MatrixXd x(as<Map<MatrixXd> >(X));

  int s = K_.length();
  std::vector<Map_MatrixXd> K;
  for(int i = 0; i < s; i++)
    K.push_back(as<Map_MatrixXd>(as<NumericMatrix>(K_[i])));

  if(verbose) Rcout << "Computing contrast matrix\n";
  HouseholderQR<MatrixXd> QR(x);
  MatrixXd Ct = QR.householderQ();
  Ct = Ct.rightCols(x.rows()-x.cols()).eval(); // aliasing !!!!

  std::vector<MatrixXd> CKCt;
  for(int i = 0; i < s; i++) {
    if(verbose) Rcout << "Computing CK[" << i << "]C'\n";
    CKCt.push_back( Ct.transpose() * K[i] * Ct );
  }

  VectorXd Cy = Ct.transpose() * y;

  double logL, logL0, normgrad;
  int niter;
  MatrixXd P(y.rows(),y.rows());
  VectorXd theta(s+1), Py(y.rows()), omega0(y.rows());
  Map_MatrixXd min_tau(as<Map_MatrixXd>(min_tau_));
  for(int i=0; i < s+1; i++) theta(i) = theta_(i);

  AIREMLn_nofix(Cy, CKCt, EMsteps, EMsteps_fail, EM_alpha, constraint, min_s2, min_tau, max_iter, eps, 
                verbose, theta, logL, logL0, niter, normgrad, P, Py, omega0, start_theta);

  List L;

  L["sigma2"] = theta(0);
  L["tau"] = theta.tail(s);
  L["logL"] = logL;
  L["logL0"] = logL0;
  L["niter"] = niter;
  L["norm_grad"] = normgrad;

  VectorXd CtPCy = Ct * Py;

  VectorXd omega(CtPCy.rows());
  omega.setZero();
  for(int i = 0; i < s; i++)
    omega.noalias() += theta(i+1)*K[i]*CtPCy;
  
  // X beta = Y - omega - sigma2 * CtPC y
  int p = x.cols();
  MatrixXd xtx( MatrixXd(p,p).setZero().selfadjointView<Lower>().rankUpdate(x.adjoint())); // X'X
  // calcul de (X'X)^{-1} X'(Y - omega)
  VectorXd beta = x.transpose() * (y - omega - theta(0)*CtPCy);
  LDLT<MatrixXd> ldlt_xtx(xtx);
  ldlt_xtx.solveInPlace(beta); 
  
  if(get_P) L["P"] = P;  
  L["Py"] = CtPCy;
  L["BLUP_omega"] = omega;
  L["BLUP_beta"] = beta;

  return L;
}

//[[Rcpp::export]]
List AIREMLn(NumericVector Y, NumericMatrix X, List K_, int EMsteps, int EMsteps_fail, double EM_alpha, bool constraint,
             double min_s2, NumericVector min_tau_, int max_iter, double eps, bool verbose,
             NumericVector theta_, bool start_theta, bool get_P) {
  Map_MatrixXd y(as<Map<MatrixXd> >(Y));
  Map_MatrixXd x(as<Map<MatrixXd> >(X));

  int s = K_.length();
  std::vector<Map_MatrixXd> K;
  for(int i = 0; i < s; i++)
    K.push_back(as<Map_MatrixXd>(as<NumericMatrix>(K_[i])));

  int p(x.cols());
  double logL, logL0, normgrad, varXbeta;
  int niter;
  VectorXd theta(s+1), Py(y.rows()), omega(y.rows()), beta(x.cols());
  MatrixXd varbeta(p,p), P(y.rows(),y.rows());
  Map_MatrixXd min_tau(as<Map_MatrixXd>(min_tau_));
  for(int i=0; i < s+1; i++) theta(i) = theta_(i);

  AIREMLn(y, x, K, EMsteps, EMsteps_fail, EM_alpha, constraint, min_s2, min_tau, max_iter, eps, verbose, theta, logL, logL0, 
            niter, normgrad, P, Py, omega, beta, varbeta, varXbeta, start_theta);

  List L;

  L["sigma2"] = theta(0);
  L["tau"] = theta.tail(s);
  L["logL"] = logL;
  L["logL0"] = logL0;
  L["niter"] = niter;
  L["norm_grad"] = normgrad;
  if(get_P) L["P"] = P;  
  L["Py"] = Py;
  L["BLUP_omega"] = omega;
  L["BLUP_beta"] = beta;
  L["varXbeta"] = varXbeta;

  return L;
}


RcppExport SEXP gg_AIREMLn_nofix(SEXP YSEXP, SEXP K_SEXP, SEXP EMstepsSEXP, SEXP EMsteps_failSEXP, SEXP EM_alphaSEXP, SEXP constraintSEXP, SEXP min_s2SEXP, SEXP min_tau_SEXP, SEXP max_iterSEXP, SEXP epsSEXP, SEXP verboseSEXP, SEXP theta_SEXP, SEXP start_thetaSEXP, SEXP get_PSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type Y(YSEXP);
    Rcpp::traits::input_parameter< List >::type K_(K_SEXP);
    Rcpp::traits::input_parameter< int >::type EMsteps(EMstepsSEXP);
    Rcpp::traits::input_parameter< int >::type EMsteps_fail(EMsteps_failSEXP);
    Rcpp::traits::input_parameter< double >::type EM_alpha(EM_alphaSEXP);
    Rcpp::traits::input_parameter< bool >::type constraint(constraintSEXP);
    Rcpp::traits::input_parameter< double >::type min_s2(min_s2SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type min_tau_(min_tau_SEXP);
    Rcpp::traits::input_parameter< int >::type max_iter(max_iterSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type theta_(theta_SEXP);
    Rcpp::traits::input_parameter< bool >::type start_theta(start_thetaSEXP);
    Rcpp::traits::input_parameter< bool >::type get_P(get_PSEXP);
    __result = Rcpp::wrap(AIREMLn_nofix(Y, K_, EMsteps, EMsteps_fail, EM_alpha, constraint, min_s2, min_tau_, max_iter, eps, verbose, theta_, start_theta, get_P));
    return __result;
END_RCPP
}

// AIREMLn_contrast
RcppExport SEXP gg_AIREMLn_contrast(SEXP YSEXP, SEXP XSEXP, SEXP K_SEXP, SEXP EMstepsSEXP, SEXP EMsteps_failSEXP, SEXP EM_alphaSEXP, SEXP constraintSEXP, SEXP min_s2SEXP, SEXP min_tau_SEXP, SEXP max_iterSEXP, SEXP epsSEXP, SEXP verboseSEXP, SEXP theta_SEXP, SEXP start_thetaSEXP, SEXP get_PSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type Y(YSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< List >::type K_(K_SEXP);
    Rcpp::traits::input_parameter< int >::type EMsteps(EMstepsSEXP);
    Rcpp::traits::input_parameter< int >::type EMsteps_fail(EMsteps_failSEXP);
    Rcpp::traits::input_parameter< double >::type EM_alpha(EM_alphaSEXP);
    Rcpp::traits::input_parameter< bool >::type constraint(constraintSEXP);
    Rcpp::traits::input_parameter< double >::type min_s2(min_s2SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type min_tau_(min_tau_SEXP);
    Rcpp::traits::input_parameter< int >::type max_iter(max_iterSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type theta_(theta_SEXP);
    Rcpp::traits::input_parameter< bool >::type start_theta(start_thetaSEXP);
    Rcpp::traits::input_parameter< bool >::type get_P(get_PSEXP);
    __result = Rcpp::wrap(AIREMLn_contrast(Y, X, K_, EMsteps, EMsteps_fail, EM_alpha, constraint, min_s2, min_tau_, max_iter, eps, verbose, theta_, start_theta, get_P));
    return __result;
END_RCPP
}

// AIREMLn
RcppExport SEXP gg_AIREMLn(SEXP YSEXP, SEXP XSEXP, SEXP K_SEXP, SEXP EMstepsSEXP, SEXP EMsteps_failSEXP, SEXP EM_alphaSEXP, SEXP constraintSEXP, SEXP min_s2SEXP, SEXP min_tau_SEXP, SEXP max_iterSEXP, SEXP epsSEXP, SEXP verboseSEXP, SEXP theta_SEXP, SEXP start_thetaSEXP, SEXP get_PSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type Y(YSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< List >::type K_(K_SEXP);
    Rcpp::traits::input_parameter< int >::type EMsteps(EMstepsSEXP);
    Rcpp::traits::input_parameter< int >::type EMsteps_fail(EMsteps_failSEXP);
    Rcpp::traits::input_parameter< double >::type EM_alpha(EM_alphaSEXP);
    Rcpp::traits::input_parameter< bool >::type constraint(constraintSEXP);
    Rcpp::traits::input_parameter< double >::type min_s2(min_s2SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type min_tau_(min_tau_SEXP);
    Rcpp::traits::input_parameter< int >::type max_iter(max_iterSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type theta_(theta_SEXP);
    Rcpp::traits::input_parameter< bool >::type start_theta(start_thetaSEXP);
    Rcpp::traits::input_parameter< bool >::type get_P(get_PSEXP);
    __result = Rcpp::wrap(AIREMLn(Y, X, K_, EMsteps, EMsteps_fail, EM_alpha, constraint, min_s2, min_tau_, max_iter, eps, verbose, theta_, start_theta, get_P));
    return __result;
END_RCPP
}
