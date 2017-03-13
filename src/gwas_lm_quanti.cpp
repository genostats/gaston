#include <Rcpp.h>
#include <RcppEigen.h>
#include "matrix4.h"
#include <ctime>
#include <cmath>

using namespace Rcpp;

// Q = une matrice avec Q'Q = Id 
// --> décomposition QR de la matrice de covariables !
//[[Rcpp::export]]
List GWAS_lm_quanti(XPtr<matrix4> pA, NumericVector mu, NumericVector Y, 
                   NumericMatrix Q_, int beg, int end) {

  Eigen::Map<Eigen::MatrixXd> Q(as<Eigen::Map<Eigen::MatrixXd> >(Q_));
  Eigen::Map<Eigen::VectorXd> y(as<Eigen::Map<Eigen::VectorXd> >(Y));

  size_t n = Y.length();
  size_t p = Q.cols();
  // résidu Y par Q
  Eigen::VectorXd z(y - Q * (Q.transpose() * y));

  NumericVector beta(end-beg+1);
  NumericVector sd_beta(end-beg+1);

  for(int i = beg; i <= end; i++) {
    Eigen::VectorXd SNP(n);
    if( std::isnan(mu(i)) || mu(i) == 0 || mu(i) == 2 ) {
      beta(i-beg) = NAN;
      sd_beta(i-beg) = NAN;
      continue;
    }

    // remplir SNP
    for(int ii = 0; ii < pA->true_ncol-1; ii++) {
      uint8_t x = pA->data[i][ii];
      for(int ss = 0; ss < 4; ss++) {
        SNP(4*ii+ss) = ((x&3) != 3)?(x&3):mu(i);
        x >>= 2;
      }
    }
    { int ii = pA->true_ncol-1;
      uint8_t x = pA->data[i][ii];
      for(int ss = 0; ss < 4 && 4*ii+ss < pA->ncol; ss++) {
        SNP(4*ii+ss) = ((x&3) != 3)?(x&3):mu(i);
        x >>= 2;
      }
    }

    // résidu
    Eigen::VectorXd g(SNP - Q * (Q.transpose() * SNP));

    // régression de z sur g
    double gg = g.squaredNorm();
    double be = g.dot(z)/gg;
    beta(i-beg) = be;
    double s2 = (z - be*g).squaredNorm()/(n-p-1); // attention au p
    sd_beta(i-beg) = sqrt(s2/gg);
  }
  List L;
  L["beta"] = beta;
  L["sd"] = sd_beta;
  return L;
}


RcppExport SEXP gg_GWAS_lm_quanti(SEXP pASEXP, SEXP muSEXP, SEXP YSEXP, SEXP Q_SEXP, SEXP begSEXP, SEXP endSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< XPtr<matrix4> >::type pA(pASEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Y(YSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Q_(Q_SEXP);
    Rcpp::traits::input_parameter< int >::type beg(begSEXP);
    Rcpp::traits::input_parameter< int >::type end(endSEXP);
    __result = Rcpp::wrap(GWAS_lm_quanti(pA, mu, Y, Q_, beg, end));
    return __result;
END_RCPP
}

