#include <Rcpp.h>
#include "diago_full.h"
#include "optimize.h"
#include "diago_wrapers.h"

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

  VectorXd Py;
  MatrixXd XViXi;
  double v;

  NumericVector res(h2.size()), s2(h2.size()), tau(h2.size());
  for(int i = 0; i < h2.size(); i++)
  {
    res(i) = diago_full_likelihood(h2(i), p, y, x, sigma, Py, v, XViXi); 
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
NumericMatrix diago_full_likelihood2(NumericVector tau, NumericVector s2, int p, NumericVector Y, NumericMatrix X, 
                                NumericVector Sigma, NumericMatrix U) {
  Map_MatrixXd y0(as<Map<MatrixXd> >(Y));
  Map_MatrixXd x0(as<Map<MatrixXd> >(X));
  Map_MatrixXd sigma(as<Map<MatrixXd> >(Sigma));
  Map_MatrixXd u(as<Map<MatrixXd> >(U));

  MatrixXd x = u.transpose() * x0;
  MatrixXd y = u.transpose() * y0;

  VectorXd Py;
  MatrixXd XViXi;

  NumericMatrix res(tau.size(), s2.size());
  for(int i = 0; i < tau.size(); i++) {
    checkUserInterrupt();
    for(int j = 0; j < s2.size(); j++)
      res(i,j) = diago_full_likelihood(tau(i), s2(j), p, y, x, sigma, Py, XViXi);
  }

  return res;
}

//[[Rcpp::export]]
List fit_full_diago(NumericVector Y, NumericMatrix X, IntegerVector p_, NumericVector Sigma, NumericMatrix U, double tol) {
  Map_MatrixXd y0(as<Map<MatrixXd> >(Y));
  Map_MatrixXd x0(as<Map<MatrixXd> >(X));
  Map_MatrixXd sigma(as<Map<MatrixXd> >(Sigma));
  Map_MatrixXd u(as<Map<MatrixXd> >(U));

  MatrixXd x = u.transpose() * x0;
  MatrixXd y = u.transpose() * y0;

  int n = sigma.rows();
  int r = x.cols();
  MatrixXd XViXi(r,r);

  List R;
  for(int i = 0; i < p_.length(); i++) {
    int p = p_(i);
    VectorXd P0y;
    double v;

    // likelihood maximization 
    par_li A;
    A.p = p; 
    A.y = &y;
    A.x = &x;
    A.sigma = &sigma;
    A.P0y = &P0y; 
    A.v = &v; 
    A.XViXi = &XViXi;

    double h2 = Brent_fmin(0, 1, wrap_li, (void *) &A, tol);
    double tau = h2*v;
    double s2 = v - tau;
    // --- end likelihood maximization

    // *********** CALCUL DES BLUPS ************************
    // Attention P0y n'est que (P0y)b, les n-p dernières composantes ! (les p premières sont nulles)
    VectorXd sigmab = sigma.bottomRows(n-p);
    VectorXd omega = h2 * sigmab.asDiagonal() * P0y;

    VectorXd z = y;
    z.tail(n-p) -= omega + (1-h2)*P0y;
    // Xb' Xb
    MatrixXd xtx( MatrixXd(r,r).setZero().selfadjointView<Lower>().rankUpdate( x.bottomRows(n-p).transpose() ));
    MatrixXd xtx0( xtx );
    MatrixXd xtxi(r,r); // et son inverse
    double d, ld;
    sym_inverse(xtx0, xtxi, d, ld, 1e-5); // détruit xtx0

    VectorXd beta(r+p);
    beta.topRows(r) = xtxi * x.bottomRows(n-p).transpose() * z.bottomRows(n-p);
    beta.bottomRows(p) = z.topRows(p) - x.topRows(p) * beta.topRows(r);
    // ************ FIN DU CALCUL DES BLUPS ******************

    // **** Calcul décomposition de la variance
    VectorXd Ut1 = u.transpose() * VectorXd::Ones(n);
    VectorXd Xbeta = x0 * beta.topRows(r) + u.leftCols(p) * beta.bottomRows(p) ;

    double psi1 = n*( p*s2 + tau*sigma.topRows(p).col(0).array().sum()  ); // n*trace(U1 Va U1') = n*trace(Va)
    psi1 -= Ut1.topRows(p).transpose() * ( tau*sigma.topRows(p).col(0).asDiagonal() )* Ut1.topRows(p);  // - 1' U1 Va U1' 1
    psi1 -= s2*Ut1.topRows(p).squaredNorm();
    psi1 /= n*(n-1);

    double psi2 = n*v*trace_of_product( xtx, XViXi ); // n*trace(U2 Xb (...) Xb' U2')
    VectorXd zz = x.bottomRows(n-p).transpose() * Ut1.bottomRows(n-p);

    psi2 -= v*zz.transpose() * XViXi * zz;
    psi2 /= n*(n-1);

    double SXbeta = Xbeta.array().sum();
    double varXbeta = (Xbeta.squaredNorm() - SXbeta*SXbeta/n)/(n-1) - psi1 - psi2;
    // **** fin décomposition !
  

    List L;
    L["sigma2"] = s2;
    L["tau"] = tau;
    L["Py"] = u.rightCols(n-p) * P0y/v;
    L["BLUP_omega"] = u.rightCols(n-p)*omega;
    L["BLUP_beta"] = beta;
    L["varbeta"] = v*XViXi;
    L["Xbeta"] = Xbeta;
    L["varXbeta"] = varXbeta;
    L["p"] = p;

    if(p_.length() > 1) R.push_back(L); else R = L;
  }
  return R;
}
