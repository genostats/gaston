#include <Rcpp.h>
#include "diago.h"
#include "optimize.h"
#include "matrix4.h"
#include "diago_wrapers.h"
#include <ctime>
#define BLOCK 20 


//[[Rcpp::export]]
List GWAS_lmm_lrt(XPtr<matrix4> pA, NumericVector mu, NumericVector Y, NumericMatrix X, 
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
  VectorXd H2(end-beg+1);
  VectorXd LRT(end-beg+1);

  // object for likelihood maximization
  par_li A;
  A.p = p; 
  A.y = &y;
  A.sigma = &sigma;
  A.P0y = &P0y; 
  A.v = &v; 
  A.XViXi = &XViXi;

  clock_t begin_t, chaviro(0);

  double h2, likelihood0;
  // on commence par le modèle null
  if(r == 1) { // pas de covariable
    begin_t = clock();
    h2 = Brent_fmin(0., 1., wrap_full_li_nc, (void *) &A, tol);
    chaviro += clock() - begin_t;
    likelihood0 = A.likelihood;
  } else {
    MatrixXd x1 = x.leftCols(r-1);
    A.x = &x1;
    begin_t = clock();
    h2 = Brent_fmin(0., 1., wrap_full_li, (void *) &A, tol);
    chaviro += clock() - begin_t;
    likelihood0 = A.likelihood;
  }

  double min_h2 = (h2 < 0.1)?0.:(h2-0.1);
  double max_h2 = (h2 > 0.9)?1.:(h2+0.1);

  // faut pas oublier d'y mettre la matrice avec une colonne de plus !!
  A.x = &x;
  for(int i = beg; i <= end; i++) {
    // remplir dernière colonne de x : récupérer SNP, multiplier par u'...
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
    x.col(r-1) = u.transpose() * SNP;

    // likelihood maximization
    begin_t = clock();
    h2 = Brent_fmin(min_h2, max_h2, wrap_full_li, (void *) &A, tol);
    chaviro += clock() - begin_t;

    if(h2 < min_h2 + tol && min_h2 > 0) { // au bas de l'intervalle
      //Rcout << "-";
      begin_t = clock();
      h2 = Brent_fmin(0., h2+tol, wrap_full_li, (void *) &A, tol);
      chaviro += clock() - begin_t;
      min_h2 = h2;
    } else if(h2 + tol > max_h2 && max_h2 < 1) { // en haut
      //Rcout << "+";
      begin_t = clock();
      h2 = Brent_fmin(h2-tol, 1., wrap_full_li, (void *) &A, tol);
      chaviro += clock() - begin_t;
      max_h2 = h2;
    }

    H2(i-beg) = h2;
    LRT(i-beg) = 2*(A.likelihood - likelihood0);

    // remettre à jour min_h2 et max_h2
    if((i - beg + 1)%BLOCK == 0) {
      min_h2 = H2.segment(i-beg-BLOCK+1,BLOCK).minCoeff();
      max_h2 = H2.segment(i-beg-BLOCK+1,BLOCK).maxCoeff();
      double delta = (max_h2 - min_h2)*0.6;
      min_h2 = (min_h2 < delta)?0.:(min_h2-delta);
      max_h2 = (max_h2 + delta > 1)?1.:(max_h2+delta);
    }

  }

  //Rcout << (float) chaviro / CLOCKS_PER_SEC << " spent in likelihood maximization\n";
  List R;
  R["h2"] = H2;
  R["LRT"] = LRT;
  return R;
}


RcppExport SEXP gg_GWAS_lmm_lrt(SEXP pASEXP, SEXP muSEXP, SEXP YSEXP, SEXP XSEXP, SEXP pSEXP, SEXP SigmaSEXP, SEXP USEXP, SEXP begSEXP, SEXP endSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< XPtr<matrix4> >::type pA(pASEXP );
        Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type Y(YSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP );
        Rcpp::traits::input_parameter< double >::type p(pSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type Sigma(SigmaSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type U(USEXP );
        Rcpp::traits::input_parameter< int >::type beg(begSEXP );
        Rcpp::traits::input_parameter< int >::type end(endSEXP );
        Rcpp::traits::input_parameter< double >::type tol(tolSEXP );
        List __result = GWAS_lmm_lrt(pA, mu, Y, X, p, Sigma, U, beg, end, tol);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}


