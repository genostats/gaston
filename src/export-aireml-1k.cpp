#include <Rcpp.h>
#include "matrix-varia.h"
#include "ai-reml-1k.h"
#include "ai-reml-1k-covar.h"

using namespace Rcpp;
typedef Map<MatrixXd> Map_MatrixXd;

// [[Rcpp::export]]
List AIREML1_nofix(NumericVector Y, NumericMatrix K_, int EMsteps, int EMsteps_fail, double EM_alpha, bool constraint, double min_s2, double min_tau, int max_iter, 
                   double eps, bool verbose, NumericVector theta0, bool start_theta, bool get_P) {
  Map_MatrixXd y(as<Map<MatrixXd> >(Y));
  Map_MatrixXd K(as<Map<MatrixXd> >(K_));

  double logL, logL0, normgrad;
  int niter;
  int n(y.rows());
  MatrixXd P(n,n);
  VectorXd Py(n), KPy(n);
  Vector2d theta;
  theta(0) = theta0(0);
  theta(1) = theta0(1);
  AIREML_nofix(y, K, EMsteps, EMsteps_fail, EM_alpha, constraint, min_s2, min_tau, max_iter, eps, 
               verbose, theta, logL, logL0, niter, normgrad, P, Py, KPy, start_theta);

  List L;
  L["sigma2"] = theta(0);
  L["tau"] = theta(1);
  L["logL"] = logL;
  L["logL0"] = logL0;
  L["niter"] = niter;
  L["norm_grad"] = normgrad;
  if(get_P) L["P"] = P;
  L["Py"] = Py; // réalise la copie 
  L["BLUP_omega"] = theta(1) * KPy;

  return L;
}

// [[Rcpp::export]]
List AIREML1_contrast(NumericVector Y, NumericMatrix X, NumericMatrix K_, int EMsteps, int EMsteps_fail, double EM_alpha, bool constraint, 
                      double min_s2, double min_tau, int max_iter, double eps, bool verbose, NumericVector theta0, bool start_theta, bool get_P) {
  Map_MatrixXd y(as<Map<MatrixXd> >(Y));
  Map_MatrixXd x(as<Map<MatrixXd> >(X));
  Map_MatrixXd K(as<Map<MatrixXd> >(K_));

  if(verbose) Rcout << "Computing contrast matrix\n";
  HouseholderQR<MatrixXd> QR(x);
  MatrixXd Ct = QR.householderQ();
  Ct = Ct.rightCols(x.rows()-x.cols()).eval(); // aliasing !!!!

  if(verbose) Rcout << "Computing CKC'\n";
  MatrixXd CKCt = Ct.transpose() * K * Ct;
  VectorXd Cy = Ct.transpose() * y;

  if(verbose) Rcout << "Fitting model\n";

  double logL, logL0, normgrad;
  int niter;
  VectorXd Py(y.rows()), KPy(y.rows());
  Vector2d theta;
  MatrixXd P(y.rows(),y.rows());
  theta(0) = theta0(0);
  theta(1) = theta0(1);
  AIREML_nofix(Cy, CKCt, EMsteps, EMsteps_fail, EM_alpha, constraint, min_s2, min_tau, max_iter, eps, 
               verbose, theta, logL, logL0, niter, normgrad, P, Py, KPy, start_theta);

  List L;
  L["sigma2"] = theta(0);
  L["tau"] = theta(1);
  L["logL"] = logL;
  L["logL0"] = logL0;
  L["niter"] = niter;
  L["norm_grad"] = normgrad;

  VectorXd CtPCy = Ct * Py;
  VectorXd omega = theta(1) * K * CtPCy;

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


// [[Rcpp::export]]
List AIREML1(NumericVector Y, NumericMatrix X, NumericMatrix K_, int EMsteps, int EMsteps_fail, double EM_alpha, bool constraint, double min_s2, 
             double min_tau, int max_iter, double eps, bool verbose, NumericVector theta0, bool start_theta, bool get_P) {
  Map_MatrixXd y(as<Map<MatrixXd> >(Y));
  Map_MatrixXd x(as<Map<MatrixXd> >(X));
  Map_MatrixXd K(as<Map<MatrixXd> >(K_));

  int p(x.cols());
  
  double logL, logL0, normgrad, varXbeta;
  int niter;
  MatrixXd P(y.rows(),y.rows());
  VectorXd Py(y.rows()), KPy(y.rows()), beta(x.cols());
  Vector2d theta;
  MatrixXd varbeta(p,p);

  theta(0) = theta0(0);
  theta(1) = theta0(1);

  AIREML1(y, x, K, EMsteps, EMsteps_fail, EM_alpha, constraint, min_s2, min_tau, max_iter, eps, verbose, theta, logL, logL0, 
           niter, normgrad, P, Py, KPy, beta, varbeta, varXbeta, start_theta);

  List L;
  L["sigma2"] = theta(0);
  L["tau"] = theta(1);
  L["logL"] = logL;
  L["logL0"] = logL0;
  L["niter"] = niter;
  L["norm_grad"] = normgrad;
  if(get_P) L["P"] = P;
  L["Py"] = Py;
  L["BLUP_omega"] = theta(1)*KPy;
  L["BLUP_beta"] = beta;
  L["varbeta"] = varbeta;
  L["varXbeta"] = varXbeta;

  return L;
}

