#include <Rcpp.h>
#include "diago_newton.h"

// h2*s + (1-h2) > 0
// h2*(s-1) > -1
// si s > 1
// h2 > 1/(1-s) donc min_h2 = max 1/(1-s) pour s > 1
// si s < 1 h2 < 1/(1-s) donc max_h2 = min 1/(1-s) pour s < 1
void min_max_h2(NumericVector Sigma, double & min_h2, double & max_h2) {
  int n = Sigma.size();
  Rcout << "n = " << n << "\n";
  max_h2 = std::numeric_limits<double>::infinity();
  min_h2 = -std::numeric_limits<double>::infinity();
  for(int i = 0; i < n; i++) {
    double s = Sigma[i];
    if(s > 1) 
      min_h2 = (min_h2 > 1/(1-s))?min_h2:1/(1-s);
    else if(s < 1)
      max_h2 = (max_h2 < 1/(1-s))?max_h2:1/(1-s);
  }
  Rcout << "min = " << min_h2 << " max = " << max_h2 << "\n";
}

//[[Rcpp::export]]
List fit_diago_newton(NumericVector Y, NumericMatrix X, IntegerVector p_, NumericVector Sigma, NumericMatrix U, double tol) {
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
   
    double v, h2 = 0.1;
    // likelihood maximization 
    double min_h2 = 0, max_h2 = 0.99;
    min_max_h2(Sigma, min_h2, max_h2);
    min_h2 += 1e-4;
    max_h2 -= 1e-4;
  Rcout << "min = " << min_h2 << " max = " << max_h2 << "\n";
    diago_likelihood_newton(h2, v, p, y, x, sigma, P0y, XViXi, min_h2, max_h2, tol, true);

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

    double s2 = (1-h2)*v, tau = h2*v;

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

RcppExport SEXP gg_fit_diago_newton(SEXP YSEXP, SEXP XSEXP, SEXP p_SEXP, SEXP SigmaSEXP, SEXP USEXP, SEXP tolSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericVector >::type Y(YSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP );
        Rcpp::traits::input_parameter< IntegerVector >::type p_(p_SEXP );
        Rcpp::traits::input_parameter< NumericVector >::type Sigma(SigmaSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type U(USEXP );
        Rcpp::traits::input_parameter< double >::type tol(tolSEXP );
        List __result = fit_diago_newton(Y, X, p_, Sigma, U, tol);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}

