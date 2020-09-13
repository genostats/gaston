#include <Rcpp.h>
#include "diago2.h"
#include "matrix4.h"
#include <cmath>

// laisser en double ça va aussi vite (plus vite ?) et ça fait vraiment
// une différence si il y a des covariables
#define scalar double
#define MATRIX MatrixXd
#define VECTOR VectorXd

//[[Rcpp::export]]
List GWAS_lmm_wald(XPtr<matrix4> pA, NumericVector mu, NumericVector Y, NumericMatrix X, 
                  int p, NumericVector Sigma, NumericMatrix U, int beg, int end, double tol) {

  int n = Sigma.size();
  int r = X.ncol();
  int max_iter = 25;

  if(Y.size() != n || X.nrow() != n || U.nrow() != n || U.ncol() != n) 
    stop("Dimensions mismatch");
    
  // conversion en float si nécessaire... sinon copie
  MATRIX y0(n, 1);
  for(int i = 0; i < n; i++)
    y0(i,0) = Y[i];

  MATRIX x0(n, r);
  for(int j = 0; j < r; j++) 
    for(int i = 0; i < n; i++)
      x0(i,j) = X(i,j);

  VECTOR sigma(n);
  for(int i = 0; i < n; i++) 
    sigma[i] = Sigma[i];

  MATRIX u(n,n);
  for(int j = 0; j < n; j++) 
    for(int i = 0; i < n; i++)
      u(i,j) = U(i,j);

  MATRIX x = u.transpose() * x0;
  MATRIX y = u.transpose() * y0;

  // Vecteur SNPs
  VECTOR SNP(n);

  // declare vectors containing result
  NumericVector H2(end-beg+1);
  NumericVector BETA(end-beg+1);
  NumericVector SDBETA(end-beg+1);

  // object for likelihood maximization
  diag_likelihood<MATRIX, VECTOR, scalar> A(p, y, x, sigma);

  scalar h2 = 0;

  for(int i = beg; i <= end; i++) {
    // if(!(i%65536)) Rcout << "i = " << i << "\n";
    if( std::isnan(mu(i)) || mu(i) == 0 || mu(i) == 2 ) {
      H2(i-beg) = NAN;
      BETA(i-beg) = NAN;
      SDBETA(i-beg) = NAN;
      continue;
    }

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
    A.X.col(r-1) = u.transpose() * SNP;

    // likelihood maximization
    h2 = (h2 > 0.9)?0.9:h2;
    A.newton_max( h2, 0, 0.99, tol, max_iter, false);
    // CALCUL DES BLUPS 
    VECTOR beta, omega;
    A.blup(h2, beta, omega, false, true);

// Rcout << "beta = " << beta.transpose() << "\n";
// Rcout << "v = " << A.v << "\n";
// Rcout << "XViX = " << A.XViX << "\n";
// Rcout << "XViX_i = " << A.XViX_i << "\n";
// Rcout << "A.d = " << A.d << "\n";
    if(A.d != 0) {
      H2(i-beg) = h2;
      BETA(i-beg) = beta(r-1);
      SDBETA(i-beg) = sqrt(A.v*A.XViX_i(r-1,r-1));
    } else {
      H2(i-beg) = NAN;
      BETA(i-beg) = NAN;
      SDBETA(i-beg) = NAN;
    }
  }

  List R;
  R["h2"] = H2;
  R["beta"] = BETA;
  R["sd"] = SDBETA;
  return R;
}


RcppExport SEXP gg_GWAS_lmm_wald(SEXP pASEXP, SEXP muSEXP, SEXP YSEXP, SEXP XSEXP, SEXP pSEXP, SEXP SigmaSEXP, SEXP USEXP, SEXP begSEXP, SEXP endSEXP, SEXP tolSEXP) {
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
        List __result = GWAS_lmm_wald(pA, mu, Y, X, p, Sigma, U, beg, end, tol);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}


