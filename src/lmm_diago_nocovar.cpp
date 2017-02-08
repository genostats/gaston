#include <Rcpp.h>
#include "diago_nocovar.h"
#include "diago_wrapers.h"
#include "optimize.h"

using namespace Rcpp;
typedef Map<MatrixXd> Map_MatrixXd;

//[[Rcpp::export]]
List diago_likelihood1_nocovar(NumericVector h2, int p, NumericVector Y, NumericVector Sigma, NumericMatrix U) {
  Map_MatrixXd y0(as<Map<MatrixXd> >(Y));
  Map_MatrixXd sigma(as<Map<MatrixXd> >(Sigma));
  Map_MatrixXd u(as<Map<MatrixXd> >(U));

  MatrixXd y = u.transpose() * y0;

  VectorXd Py;
  double v;

  NumericVector res(h2.size()), s2(h2.size()), tau(h2.size());
  for(int i = 0; i < h2.size(); i++) {
    res(i) = diago_likelihood(h2(i), p, y, sigma, Py, v); 
    tau(i) = h2(i)*v;
    s2(i) = (1-h2(i))*v;
  }

  List L;
  L["tau"] = tau;
  L["sigma2"] = s2;
  L["likelihood"] = res;
  return L;
}

//[[Rcpp::export]]
NumericMatrix diago_likelihood2_nocovar(NumericVector tau, NumericVector s2, int p, NumericVector Y, 
                                NumericVector Sigma, NumericMatrix U) {
  Map_MatrixXd y0(as<Map<MatrixXd> >(Y));
  Map_MatrixXd sigma(as<Map<MatrixXd> >(Sigma));
  Map_MatrixXd u(as<Map<MatrixXd> >(U));

  MatrixXd y = u.transpose() * y0;

  VectorXd Py;

  NumericMatrix res(tau.size(), s2.size());
  for(int i = 0; i < tau.size(); i++) {
    checkUserInterrupt();
    for(int j = 0; j < s2.size(); j++)
      res(i,j) = diago_likelihood(tau(i), s2(j), p, y, sigma, Py);
  }

  return res;
}

//[[Rcpp::export]]
List fit_diago_nocovar(NumericVector Y, IntegerVector p_, NumericVector Sigma, NumericMatrix U, double tol) {
  Map_MatrixXd y0(as<Map<MatrixXd> >(Y));
  Map_MatrixXd sigma(as<Map<MatrixXd> >(Sigma));
  Map_MatrixXd u(as<Map<MatrixXd> >(U));

  MatrixXd y = u.transpose() * y0;

  int n = sigma.rows();

  List R;
  for(int i = 0; i < p_.length(); i++) {
    int p = p_(i);
    VectorXd P0y;
    double v;

    // likelihood maximization 
    par_li A;
    A.p = p; 
    A.y = &y;
    A.sigma = &sigma;
    A.P0y = &P0y; 
    A.v = &v; 

    double h2 = Brent_fmin(0, 1, wrap_li_nc, (void *) &A, tol);
    double tau = h2*v;
    double s2 = v - tau;
    // --- end likelihood maximization

    // *********** CALCUL DES BLUPS ************************
    // Attention P0y n'est que (P0y)b, les n-p dernières composantes ! (les p premières sont nulles)
    VectorXd sigmab = sigma.bottomRows(n-p);
    VectorXd omega = h2 * sigmab.asDiagonal() * P0y;

    VectorXd beta = y.topRows(p);
    // ************ FIN DU CALCUL DES BLUPS ******************

    // **** Calcul décomposition de la variance
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
    L["Py"] = u.rightCols(n-p) * P0y/v;
    L["BLUP_omega"] = u.rightCols(n-p)*omega;
    L["BLUP_beta"] = beta;
    L["Xbeta"] = Xbeta;
    L["varXbeta"] = varXbeta;
    L["p"] = p;

    if(p_.length() > 1) R.push_back(L); else R = L;
  }
  return R;
}

RcppExport SEXP gg_diago_likelihood1_nocovar(SEXP h2SEXP, SEXP pSEXP, SEXP YSEXP, SEXP SigmaSEXP, SEXP USEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericVector >::type h2(h2SEXP );
        Rcpp::traits::input_parameter< int >::type p(pSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type Y(YSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type Sigma(SigmaSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type U(USEXP );
        List __result = diago_likelihood1_nocovar(h2, p, Y, Sigma, U);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}

RcppExport SEXP gg_diago_likelihood2_nocovar(SEXP tauSEXP, SEXP s2SEXP, SEXP pSEXP, SEXP YSEXP, SEXP SigmaSEXP, SEXP USEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericVector >::type tau(tauSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type s2(s2SEXP );
        Rcpp::traits::input_parameter< int >::type p(pSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type Y(YSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type Sigma(SigmaSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type U(USEXP );
        NumericMatrix __result = diago_likelihood2_nocovar(tau, s2, p, Y, Sigma, U);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}

RcppExport SEXP gg_fit_diago_nocovar(SEXP YSEXP, SEXP p_SEXP, SEXP SigmaSEXP, SEXP USEXP, SEXP tolSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericVector >::type Y(YSEXP );
        Rcpp::traits::input_parameter< IntegerVector >::type p_(p_SEXP );
        Rcpp::traits::input_parameter< NumericVector >::type Sigma(SigmaSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type U(USEXP );
        Rcpp::traits::input_parameter< double >::type tol(tolSEXP );
        List __result = fit_diago_nocovar(Y, p_, Sigma, U, tol);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}


