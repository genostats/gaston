#include <Rcpp.h>
#include "diago.h"
#include <RcppParallel.h>
#include "optimize.h"
#include "diago_wrapers.h"
#include "matrix4.h"
#include <ctime>

#define BLOCK 20 

using namespace Rcpp;
using namespace RcppParallel;

struct paraLRT : public Worker {
  // input
  uint8_t ** const data;
  size_t ncol;
  size_t true_ncol;
  size_t firstsnp;
  const double * mu;
  const MatrixXd x_, y;
  const Map_MatrixXd sigma, u;

  int n, r, p;
  double tol;
  double min_h2_;
  double max_h2_;  

  // output
  double * H2, * LIK;

  paraLRT(uint8_t ** const data, size_t ncol, size_t true_ncol, size_t firstsnp, const double * mu, const MatrixXd & x, const MatrixXd & y, const Map_MatrixXd & sigma, 
           const Map_MatrixXd & u, int n, int r, int p, double tol, double min_h2, double max_h2, double * H2, double * LIK) :
           data(data), ncol(ncol), true_ncol(true_ncol), firstsnp(firstsnp), x_(x), y(y), sigma(sigma), 
           u(u), n(n), r(r), p(p), tol(tol), min_h2_(min_h2), max_h2_(max_h2), H2(H2), LIK(LIK) {}

  void operator()(size_t beg, size_t end) {
    // declare vectors used in the loop
    VectorXd P0y(n-p);
    VectorXd beta(r);
    VectorXd sigmab(n-p);
    VectorXd omega(n-p);
    VectorXd z(n); // possible de se contenter de n-p et d'éviter les tail
    MatrixXd XViXi(r,r);
    double v;

    // recopier x !!
    MatrixXd x(n,r);
    x.leftCols(r-1) = x_.leftCols(r-1);
    // et min max h2 qui sont updatés dans chaque morceau a priori...
    double min_h2 = min_h2_, max_h2 = max_h2_;

    par_li A;
    A.p = p; 
    A.y = &y;
    A.x = &x;
    A.sigma = &sigma;
    A.P0y = &P0y; 
    A.v = &v; 
    A.XViXi = &XViXi;

    // dernière colonne de x
    VectorXd SNP(n);

    for(int i = beg; i < end; i++) {
      // remplir dernière colonne de x : récupérer SNP, multiplier par u'...
      for(int ii = 0; ii < true_ncol-1; ii++) {
        uint8_t xx = data[firstsnp+i][ii];
        for(int ss = 0; ss < 4; ss++) {
          SNP(4*ii+ss) = ((xx&3) != 3)?(xx&3):mu[i];
          xx >>= 2;
        }
      }
      { int ii = true_ncol-1;
        uint8_t xx = data[firstsnp+i][ii];
        for(int ss = 0; ss < 4 && 4*ii+ss < ncol; ss++) {
          SNP(4*ii+ss) = ((xx&3) != 3)?(xx&3):mu[i];
          xx >>= 2;
        }
      }
      x.col(r-1) = u.transpose() * SNP;
      // likelihood maximization
      double h2 = Brent_fmin(min_h2, max_h2, wrap_full_li, (void *) &A, tol);
      
      if(h2 < min_h2 + tol && min_h2 > 0) { // au bas de l'intervalle
        //std::cout << "-";
        h2 = Brent_fmin(0., h2+tol, wrap_full_li, (void *) &A, tol);
        min_h2 = h2;
      } else if(h2 + tol > max_h2 && max_h2 < 1) { // en haut
        //std::cout << "+";
        h2 = Brent_fmin(h2-tol, 1., wrap_full_li, (void *) &A, tol);
        max_h2 = h2;
      }
      H2[i] = h2;
      LIK[i] = A.likelihood;

      // remettre à jour min_h2 et max_h2
      if((i - beg + 1)%BLOCK == 0) {
        for(int ii = 0; ii < BLOCK; ii++) {
          min_h2 = (H2[i-ii] < min_h2)?H2[i-ii]:min_h2;
          max_h2 = (H2[i-ii] > max_h2)?H2[i-ii]:max_h2;
        }
        double delta = (max_h2 - min_h2)*0.6;
        min_h2 = (min_h2 < delta)?0.:(min_h2-delta);
        max_h2 = (max_h2 + delta > 1)?1.:(max_h2+delta);
      }
    }
  }

};


//[[Rcpp::export]]
List GWAS_lmm_lrt_mt(XPtr<matrix4> pA, NumericVector mu, NumericVector Y, NumericMatrix X, 
                  int p, NumericVector Sigma, NumericMatrix U, int beg, int end, double tol) {
  Map_MatrixXd y0(as<Map<MatrixXd> >(Y));
  Map_MatrixXd x0(as<Map<MatrixXd> >(X));
  Map_MatrixXd sigma(as<Map<MatrixXd> >(Sigma));
  Map_MatrixXd u(as<Map<MatrixXd> >(U));

  MatrixXd x = u.transpose() * x0;
  MatrixXd y = u.transpose() * y0;
 
  int n = sigma.rows();
  int r = x.cols();
  MatrixXd XViXi(r,r);

  double v;

  // declare vectors used in the loop
  VectorXd P0y(n-p);
  VectorXd beta(r);
  VectorXd sigmab(n-p);
  VectorXd omega(n-p);
  VectorXd z(n); // possible de se contenter de n-p et d'éviter les tail

  // dernière colonne de r
  VectorXd SNP(n);

  // declare vectors containing result
  NumericVector H2(end-beg+1);
  NumericVector LIK(end-beg+1);

  // object for likelihood maximization
  par_li A;
  A.p = p;
  A.y = &y;
  A.sigma = &sigma;
  A.P0y = &P0y;
  A.v = &v;
  A.XViXi = &XViXi;

  paraLRT XX(pA->data, pA->ncol, pA->true_ncol, beg, mu.begin(), x, y, sigma, u, n, r, p, tol, 0., 1.,
         H2.begin(), LIK.begin());

  double h2, likelihood0;
  // on commence par le modèle null
  if(r == 1) { // pas de covariable
    h2 = Brent_fmin(0., 1., wrap_full_li_nc, (void *) &A, tol);
    likelihood0 = A.likelihood;
  } else {
    MatrixXd x1 = x.leftCols(r-1);
    A.x = &x1;
    h2 = Brent_fmin(0., 1., wrap_full_li, (void *) &A, tol);
    likelihood0 = A.likelihood;
  }
  XX.min_h2_ = (h2 < 0.1)?0.:(h2-0.1);
  XX.max_h2_ = (h2 > 0.9)?1.:(h2+0.1);

  int lea = (end-beg+1)/100; 
  if(lea < 20) {
    parallelFor(0, end-beg+1, XX, end-beg+1); // on ne parallélise pas un truc si petit
  } else {
    parallelFor(0, lea, XX, lea); // on fait le début pour estimer min et max h2
    double min_h2 = 1, max_h2 = 0;
    for(int i = 0; i < lea; i++) {
      min_h2 = (H2(i) < min_h2)?H2(i):min_h2;
      max_h2 = (H2(i) > max_h2)?H2(i):max_h2;
    } 
    XX.min_h2_ = min_h2; // et on les utilise pour le reste
    XX.max_h2_ = max_h2;
    parallelFor(lea, end-beg+1, XX);
  }

  List R;
  R["h2"] = H2;
  R["LRT"] = 2*(LIK-likelihood0);
  return R;
}

RcppExport SEXP gg_GWAS_lmm_lrt_mt(SEXP pASEXP, SEXP muSEXP, SEXP YSEXP, SEXP XSEXP, SEXP pSEXP, SEXP SigmaSEXP, SEXP USEXP, SEXP begSEXP, SEXP endSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< XPtr<matrix4> >::type pA(pASEXP );
        Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type Y(YSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP );
        Rcpp::traits::input_parameter< int >::type p(pSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type Sigma(SigmaSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type U(USEXP );
        Rcpp::traits::input_parameter< int >::type beg(begSEXP );
        Rcpp::traits::input_parameter< int >::type end(endSEXP );
        Rcpp::traits::input_parameter< double >::type tol(tolSEXP );
        List __result = GWAS_lmm_lrt_mt(pA, mu, Y, X, p, Sigma, U, beg, end, tol);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}

