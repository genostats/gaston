#include <Rcpp.h>
#include "diago2_full.h"
#include "diago2_full_nocovar.h"
#include "matrix4.h"

// laisser en double ça va aussi vite (plus vite ?) et ça fait vraiment
// une différence si il y a des covariables
#define scalar double
#define MATRIX MatrixXd
#define VECTOR VectorXd

//[[Rcpp::export]]
List GWAS_lmm_lrt(XPtr<matrix4> pA, NumericVector mu, NumericVector Y, NumericMatrix X, 
                  int p, NumericVector Sigma, NumericMatrix U, int beg, int end, double tol) {

  int n = Sigma.size();
  int r = X.ncol();
  int max_iter = 25;

  if(Y.size() != n || X.nrow() != n || U.nrow() != n || U.ncol() != n)
    stop("Dimensions mismatch");

  // conversion en float...
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

  // Zecteur SNPs
  VECTOR SNP(n);

  // declare vectors containing result
  NumericVector H2(end-beg+1);
  NumericVector LRT(end-beg+1);

  scalar h2 = 0;
  scalar likelihood0;

  // on commence par le modèle null
  if(r == 1) { 
    // object for likelihood maximization
    diag_full_likelihood_nocovar<MATRIX, VECTOR, scalar> A(p, y, sigma);
    A.newton_max(h2, 0, 0.99, tol, max_iter, false);
    likelihood0 = A.likelihood();
  } else { 
    MATRIX x1 = x.leftCols(r-1);
    // object for likelihood maximization
    diag_full_likelihood<MATRIX, VECTOR, scalar> A(p, y, x1, sigma);
    A.newton_max(h2, 0, 0.99, tol, max_iter, false);
    likelihood0 = A.likelihood();
  }
  // object for likelihood maximization
  diag_full_likelihood<MATRIX, VECTOR, scalar> A(p, y, x, sigma);

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
    A.X.col(r-1) = u.transpose() * SNP;

    // likelihood maximization
    h2 = (h2 > 0.9)?0.9:h2;
    A.newton_max( h2, 0, 0.99, tol, max_iter, false);

    H2(i-beg) = h2;
    LRT(i-beg) = 2*(A.likelihood() - likelihood0);
  }

  List R;
  R["h2"] = H2;
  R["LRT"] = LRT;
  return R;
}

