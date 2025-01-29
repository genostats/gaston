#include <Rcpp.h>
#include "matrix4.h"
#include "logit_model.h"
#include <ctime>
#include <cmath>
#include <iostream>
#define BLOCK 20 

//[[Rcpp::export]]
List GWAS_logit_wald_f(XPtr<matrix4> pA, NumericVector mu, NumericVector Y, NumericMatrix X, 
                       int beg, int end, double tol) {
  int n = Y.size();
  int r = X.ncol();

  // recopiage des matrices... en float
  MatrixXf y(n,1);
  MatrixXf x(n,r);
  for(int i = 0; i < n; i++) y(i,0) = (float) Y[i];

  for(int i = 0; i < n; i++) 
    for(int j = 0; j < r; j++)
      x(i,j) = (float) X(i,j);

  // declare vectors containing result
  NumericVector BETA(end-beg+1);
  NumericVector SDBETA(end-beg+1);

  VectorXf beta(r);
  beta.setZero();
  for(int i = beg; i <= end; i++) { 
    if( std::isnan(mu(i)) || mu(i) == 0 || mu(i) == 2 ) {
      BETA(i-beg) = NAN;
      SDBETA(i-beg) = NAN;
      continue;
    }
    // remplir dernière colonne de x par génotype au SNP (manquant -> mu)
    for(int ii = 0; ii < pA->true_ncol-1; ii++) {
      uint8_t xx = pA->data[i][ii];
      for(int ss = 0; ss < 4; ss++) {
        x(4*ii+ss, r-1) = ((xx&3) != 3)?(xx&3):mu(i);
        xx >>= 2;
      }
    }
    { int ii = pA->true_ncol-1;
      uint8_t xx = pA->data[i][ii];
      for(int ss = 0; ss < 4 && 4*ii+ss < pA->ncol; ss++) {
        x(4*ii+ss, r-1) = ((xx&3) != 3)?(xx&3):mu(i);
        xx >>= 2;
      }
    }

    MatrixXf varbeta(r,r);
    logistic_model2<float>(y, x, beta, varbeta, tol);

    BETA(i-beg) = beta(r-1);
    SDBETA(i-beg) = sqrt(varbeta(r-1,r-1));
  }

  //cout << (float) chaviro / CLOCKS_PER_SEC << "\n";
  List R;
  R["beta"] = BETA;
  R["sd"] = SDBETA;
  return R;
}

