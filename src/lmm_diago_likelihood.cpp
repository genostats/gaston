#include <Rcpp.h>
#include "diago2.h"
#include "diago2_nocovar.h"

using namespace Rcpp;
typedef Map<MatrixXd> Map_MatrixXd;

//[[Rcpp::export]]
List diago_likelihood1(NumericVector h2, int p, NumericVector Y, NumericMatrix X, NumericVector Sigma, NumericMatrix U) {
  Map_MatrixXd y0(as<Map<MatrixXd> >(Y));
  Map_MatrixXd x0(as<Map<MatrixXd> >(X));
  Map_MatrixXd sigma(as<Map<MatrixXd> >(Sigma));
  Map_MatrixXd u(as<Map<MatrixXd> >(U));

  MatrixXd x = u.transpose() * x0;
  MatrixXd y = u.transpose() * y0;

  diag_likelihood<MatrixXd, VectorXd, double> A(p, y, x, sigma);

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
NumericMatrix diago_likelihood2(NumericVector tau, NumericVector s2, int p, NumericVector Y, NumericMatrix X, 
                                NumericVector Sigma, NumericMatrix U) {
  Map_MatrixXd y0(as<Map<MatrixXd> >(Y));
  Map_MatrixXd x0(as<Map<MatrixXd> >(X));
  Map_MatrixXd sigma(as<Map<MatrixXd> >(Sigma));
  Map_MatrixXd u(as<Map<MatrixXd> >(U));

  MatrixXd x = u.transpose() * x0;
  MatrixXd y = u.transpose() * y0;

  diag_likelihood<MatrixXd, VectorXd, double> A(p, y, x, sigma);

  NumericMatrix res(tau.size(), s2.size());
  for(int i = 0; i < tau.size(); i++) {
    checkUserInterrupt();
    for(int j = 0; j < s2.size(); j++)
      res(i,j) = A.likelihood(tau(i), s2(j));
  }

  return res;
}

//[[Rcpp::export]]
List diago_likelihood1_nocovar(NumericVector h2, int p, NumericVector Y, NumericVector Sigma, NumericMatrix U) {
  Map_MatrixXd y0(as<Map<MatrixXd> >(Y));
  Map_MatrixXd sigma(as<Map<MatrixXd> >(Sigma));
  Map_MatrixXd u(as<Map<MatrixXd> >(U));

  MatrixXd y = u.transpose() * y0;

  diag_likelihood_nocovar<MatrixXd, VectorXd, double> A(p, y, sigma);

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
NumericMatrix diago_likelihood2_nocovar(NumericVector tau, NumericVector s2, int p, NumericVector Y, 
                                NumericVector Sigma, NumericMatrix U) {
  Map_MatrixXd y0(as<Map<MatrixXd> >(Y));
  Map_MatrixXd sigma(as<Map<MatrixXd> >(Sigma));
  Map_MatrixXd u(as<Map<MatrixXd> >(U));

  MatrixXd y = u.transpose() * y0;

  diag_likelihood_nocovar<MatrixXd, VectorXd, double> A(p, y, sigma);

  NumericMatrix res(tau.size(), s2.size());
  for(int i = 0; i < tau.size(); i++) {
    checkUserInterrupt();
    for(int j = 0; j < s2.size(); j++)
      res(i,j) = A.likelihood(tau(i), s2(j));
  }

  return res;
}

