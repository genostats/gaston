#include <Rcpp.h>
#include "diago2.h"
#include "diago2_nocovar.h"
#include "lmm_diago_min_max_h2.h"

using namespace Rcpp;

List fit_diago(NumericVector Y, NumericMatrix X, IntegerVector p_, NumericVector Sigma, NumericMatrix U, double min_h2, double max_h2, double tol, double verbose, bool brent) {
  Map_MatrixXd y0(as<Map<MatrixXd> >(Y));
  Map_MatrixXd x0(as<Map<MatrixXd> >(X));
  Map<VectorXd> sigma(as<Map<VectorXd> >(Sigma));
  Map_MatrixXd u(as<Map<MatrixXd> >(U));

  MatrixXd x = u.transpose() * x0;
  MatrixXd y = u.transpose() * y0;

  int n = sigma.rows();
  int r = x.cols();  
  int max_iter = 10;

  List R;

  min_max_h2(Sigma, min_h2, max_h2);
  if(verbose) Rcout << "Optimization in interval [" << min_h2 << ", " << max_h2 << "]" << std::endl;

  for(int i = 0; i < p_.length(); i++) {
    int p = p_(i);
    if(verbose) Rcout << "Optimizing with p = " << p << "\n";
   
    // likelihood maximization 
    double h2 = min_h2;
    diag_likelihood<MatrixXd, VectorXd, double> A(p, y, x, sigma);
    if(brent) 
      h2 = A.Brent_fmax(min_h2, max_h2, tol);
    else
      A.newton_max(h2, min_h2, max_h2, tol, max_iter, verbose);

    // calcul blups transféré dans la classe diag_likelihood !...
    VectorXd beta, omega;
    A.blup(h2, beta, omega, false, false);
    double s2 = (1-h2)*A.v, tau = h2*A.v;

    // **** Calcul décomposition de la variance gardé ici (on a besoin de la matrice u)
    VectorXd Ut1 = u.transpose() * VectorXd::Ones(n);
    VectorXd Xbeta = x0 * beta.topRows(r) + u.leftCols(p) * beta.bottomRows(p) ;

    double psi1 = n*( p*s2 + tau*sigma.topRows(p).col(0).array().sum()  ); // n*trace(U1 Va U1') = n*trace(Va)
    psi1 -= Ut1.topRows(p).transpose() * ( tau*sigma.topRows(p).col(0).asDiagonal() )* Ut1.topRows(p);  // - 1' U1 Va U1' 1
    psi1 -= s2*Ut1.topRows(p).squaredNorm();
    psi1 /= n*(n-1);

    double psi2 = n*A.v*trace_of_product( A.xtx, A.XViX_i ); // n*trace(U2 Xb (...) Xb' U2')
    VectorXd zz = x.bottomRows(n-p).transpose() * Ut1.bottomRows(n-p);

    psi2 -= A.v*zz.transpose() * A.XViX_i * zz;
    psi2 /= n*(n-1);

    double SXbeta = Xbeta.array().sum();
    double varXbeta = (Xbeta.squaredNorm() - SXbeta*SXbeta/n)/(n-1) - psi1 - psi2;
    // **** fin décomposition !
  

    List L;
    L["sigma2"] = s2;
    L["tau"] = tau;
    L["Py"] = u.rightCols(n-p) * A.P0y/A.v;
    L["BLUP_omega"] = u.rightCols(n-p)*omega;
    L["BLUP_beta"] = beta;
    L["varbeta"] = A.v*A.XViX_i;
    L["Xbeta"] = Xbeta;
    L["varXbeta"] = varXbeta;
    L["p"] = p;

    if(p_.length() > 1) R.push_back(L); else R = L;
  }
  return R;
}

//[[Rcpp::export]]
List fit_diago_nocovar(NumericVector Y, IntegerVector p_, NumericVector Sigma, NumericMatrix U, double min_h2, double max_h2, double tol, double verbose, double brent) {
  Map_MatrixXd y0(as<Map<MatrixXd> >(Y));
  Map<VectorXd> sigma(as<Map<VectorXd> >(Sigma));
  Map_MatrixXd u(as<Map<MatrixXd> >(U));

  MatrixXd y = u.transpose() * y0;

  int n = sigma.rows();
  int max_iter = 10;
  List R;

  min_max_h2(Sigma, min_h2, max_h2);
  if(verbose) Rcout << "Optimization in interval [" << min_h2 << ", " << max_h2 << "]" << std::endl;

  for(int i = 0; i < p_.length(); i++) {
    int p = p_(i);
    if(verbose) Rcout << "Optimizing with p = " << p << "\n";
   
    // likelihood maximization 
    double h2 = min_h2;
    diag_likelihood_nocovar<MatrixXd, VectorXd, double> A(p, y, sigma);
    if(brent) 
      h2 = A.Brent_fmax(min_h2, max_h2, tol);
    else
      A.newton_max(h2, min_h2, max_h2, tol, max_iter, verbose);

    // calcul blups transféré dans la classe diag_likelihood !...
    VectorXd beta, omega;
    A.blup(h2, beta, omega, false);
    double s2 = (1-h2)*A.v, tau = h2*A.v;

    // **** Calcul décomposition de la variance gardé ici (on a besoin de la matrice u)
    VectorXd Ut1 = u.transpose() * VectorXd::Ones(n);
    VectorXd Xbeta = u.leftCols(p) * beta.bottomRows(p) ;

    double psi1 = n*( p*s2 + tau*sigma.topRows(p).col(0).array().sum()  ); // n*trace(U1 Va U1') = n*trace(Va)
    psi1 -= Ut1.topRows(p).transpose() * ( tau*sigma.topRows(p).col(0).asDiagonal() )* Ut1.topRows(p);  // - 1' U1 Va U1' 1
    psi1 -= s2*Ut1.topRows(p).squaredNorm();
    psi1 /= n*(n-1);

    double SXbeta = Xbeta.array().sum();
    double varXbeta = (Xbeta.squaredNorm() - SXbeta*SXbeta/n)/(n-1) - psi1;
    // **** fin décomposition !
  

    List L;
    L["sigma2"] = s2;
    L["tau"] = tau;
    L["Py"] = u.rightCols(n-p) * A.P0y/A.v;
    L["BLUP_omega"] = u.rightCols(n-p)*omega;
    L["BLUP_beta"] = beta;
    L["Xbeta"] = Xbeta;
    L["varXbeta"] = varXbeta;
    L["p"] = p;

    if(p_.length() > 1) R.push_back(L); else R = L;
  }
  return R;
}

RcppExport SEXP gg_fit_diago(SEXP YSEXP, SEXP XSEXP, SEXP p_SEXP, SEXP SigmaSEXP, SEXP USEXP, SEXP min_h2_SEXP, SEXP max_h2_SEXP, SEXP tolSEXP, SEXP verbose_SEXP, SEXP brent_SEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericVector >::type Y(YSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP );
        Rcpp::traits::input_parameter< IntegerVector >::type p_(p_SEXP );
        Rcpp::traits::input_parameter< NumericVector >::type Sigma(SigmaSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type U(USEXP );
        Rcpp::traits::input_parameter< double >::type min_h2(min_h2_SEXP );
        Rcpp::traits::input_parameter< double >::type max_h2(max_h2_SEXP );
        Rcpp::traits::input_parameter< double >::type tol(tolSEXP );
        Rcpp::traits::input_parameter< bool >::type verbose(verbose_SEXP );
        Rcpp::traits::input_parameter< bool >::type brent(brent_SEXP );
        List __result = fit_diago(Y, X, p_, Sigma, U, min_h2, max_h2, tol, verbose, brent);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}

RcppExport SEXP gg_fit_diago_nocovar(SEXP YSEXP, SEXP p_SEXP, SEXP SigmaSEXP, SEXP USEXP, SEXP min_h2_SEXP, SEXP max_h2_SEXP, SEXP tolSEXP, SEXP verbose_SEXP, SEXP brent_SEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericVector >::type Y(YSEXP );
        Rcpp::traits::input_parameter< IntegerVector >::type p_(p_SEXP );
        Rcpp::traits::input_parameter< NumericVector >::type Sigma(SigmaSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type U(USEXP );
        Rcpp::traits::input_parameter< double >::type min_h2(min_h2_SEXP );
        Rcpp::traits::input_parameter< double >::type max_h2(max_h2_SEXP );
        Rcpp::traits::input_parameter< double >::type tol(tolSEXP );
        Rcpp::traits::input_parameter< bool >::type verbose(verbose_SEXP );
        Rcpp::traits::input_parameter< bool >::type brent(brent_SEXP );
        List __result = fit_diago_nocovar(Y, p_, Sigma, U, min_h2, max_h2, tol, verbose, brent);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}

