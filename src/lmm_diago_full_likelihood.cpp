#include <Rcpp.h>
#include "diago2_full.h"
#include "diago2_full_nocovar.h"

using namespace Rcpp;
typedef Map<MatrixXd> Map_MatrixXd;

//[[Rcpp::export]]
List diago_full_likelihood1(NumericVector h2, int p, NumericVector Y, NumericMatrix X, NumericVector Sigma, NumericMatrix U) {
  Map_MatrixXd y0(as<Map<MatrixXd> >(Y));
  Map_MatrixXd x0(as<Map<MatrixXd> >(X));
  Map_MatrixXd sigma(as<Map<MatrixXd> >(Sigma));
  Map_MatrixXd u(as<Map<MatrixXd> >(U));

  MatrixXd x = u.transpose() * x0;
  MatrixXd y = u.transpose() * y0;

  diag_full_likelihood<MatrixXd, VectorXd, double> A(p, y, x, sigma);

  NumericVector res(h2.size()), s2(h2.size()), tau(h2.size());
  for(int i = 0; i < h2.size(); i++) {
    res(i) = A.f(h2(i)); 
    tau(i) = h2(i)*A.v;
    s2(i) = (1-h2(i))*A.v;
  }

  List L;
  L["tau"] = tau;
  L["sigma2"] = s2;
  L["likelihood"] = res;
  return L;
}

//[[Rcpp::export]]
NumericMatrix diago_full_likelihood2(NumericVector tau, NumericVector s2, int p, NumericVector Y, NumericMatrix X, 
                                NumericVector Sigma, NumericMatrix U) {
  Map_MatrixXd y0(as<Map<MatrixXd> >(Y));
  Map_MatrixXd x0(as<Map<MatrixXd> >(X));
  Map_MatrixXd sigma(as<Map<MatrixXd> >(Sigma));
  Map_MatrixXd u(as<Map<MatrixXd> >(U));

  MatrixXd x = u.transpose() * x0;
  MatrixXd y = u.transpose() * y0;

  diag_full_likelihood<MatrixXd, VectorXd, double> A(p, y, x, sigma);

  NumericMatrix res(tau.size(), s2.size());
  for(int i = 0; i < tau.size(); i++) {
    checkUserInterrupt();
    for(int j = 0; j < s2.size(); j++)
      res(i,j) = A.likelihood(tau(i), s2(j));
  }

  return res;
}

//[[Rcpp::export]]
List diago_full_likelihood1_nocovar(NumericVector h2, int p, NumericVector Y, NumericVector Sigma, NumericMatrix U) {
  Map_MatrixXd y0(as<Map<MatrixXd> >(Y));
  Map_MatrixXd sigma(as<Map<MatrixXd> >(Sigma));
  Map_MatrixXd u(as<Map<MatrixXd> >(U));

  MatrixXd y = u.transpose() * y0;

  diag_full_likelihood_nocovar<MatrixXd, VectorXd, double> A(p, y, sigma);

  NumericVector res(h2.size()), s2(h2.size()), tau(h2.size());
  for(int i = 0; i < h2.size(); i++) {
    res(i) = A.f(h2(i));
    tau(i) = h2(i)*A.v;
    s2(i) = (1-h2(i))*A.v;
  }

  List L;
  L["tau"] = tau;
  L["sigma2"] = s2;
  L["likelihood"] = res;
  return L;
}

//[[Rcpp::export]]
NumericMatrix diago_full_likelihood2_nocovar(NumericVector tau, NumericVector s2, int p, NumericVector Y, 
                                NumericVector Sigma, NumericMatrix U) {
  Map_MatrixXd y0(as<Map<MatrixXd> >(Y));
  Map_MatrixXd sigma(as<Map<MatrixXd> >(Sigma));
  Map_MatrixXd u(as<Map<MatrixXd> >(U));

  MatrixXd y = u.transpose() * y0;

  diag_full_likelihood_nocovar<MatrixXd, VectorXd, double> A(p, y, sigma);

  NumericMatrix res(tau.size(), s2.size());
  for(int i = 0; i < tau.size(); i++) {
    checkUserInterrupt();
    for(int j = 0; j < s2.size(); j++)
      res(i,j) = A.likelihood(tau(i), s2(j));
  }

  return res;
}


RcppExport SEXP gg_diago_full_likelihood1_nocovar(SEXP h2SEXP, SEXP pSEXP, SEXP YSEXP, SEXP SigmaSEXP, SEXP USEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericVector >::type h2(h2SEXP );
        Rcpp::traits::input_parameter< int >::type p(pSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type Y(YSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type Sigma(SigmaSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type U(USEXP );
        List __result = diago_full_likelihood1_nocovar(h2, p, Y, Sigma, U);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}

RcppExport SEXP gg_diago_full_likelihood2_nocovar(SEXP tauSEXP, SEXP s2SEXP, SEXP pSEXP, SEXP YSEXP, SEXP SigmaSEXP, SEXP USEXP) {
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
        NumericMatrix __result = diago_full_likelihood2_nocovar(tau, s2, p, Y, Sigma, U);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}


RcppExport SEXP gg_diago_full_likelihood1(SEXP h2SEXP, SEXP pSEXP, SEXP YSEXP, SEXP XSEXP, SEXP SigmaSEXP, SEXP USEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericVector >::type h2(h2SEXP );
        Rcpp::traits::input_parameter< int >::type p(pSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type Y(YSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type Sigma(SigmaSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type U(USEXP );
        List __result = diago_full_likelihood1(h2, p, Y, X, Sigma, U);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}

RcppExport SEXP gg_diago_full_likelihood2(SEXP tauSEXP, SEXP s2SEXP, SEXP pSEXP, SEXP YSEXP, SEXP XSEXP, SEXP SigmaSEXP, SEXP USEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericVector >::type tau(tauSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type s2(s2SEXP );
        Rcpp::traits::input_parameter< int >::type p(pSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type Y(YSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type Sigma(SigmaSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type U(USEXP );
        NumericMatrix __result = diago_full_likelihood2(tau, s2, p, Y, X, Sigma, U);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}


