#include <iostream>
#include "loubar.h"
using namespace Rcpp;

/******** produit lou x bar ****/

void loubar(const lou & A, const bar & X, bar & R) {
  if(A.ncol != X.n || A.nrow != R.n) {
    Rcerr << "dim mismatch in loubar\n";
    return;
  }
  R.zeros();
  size_t off = 0;
  for(size_t j = 0; j < A.ncol; j++) {
    double Xj = X.data[j];
    for(size_t i = 0; i < A.nrow; i++) 
      R.data[i] += A.data[off++]*Xj;
  }
}


bar loubar(const lou & A, const bar & X) { 
  bar R(A.nrow); 
  loubar(A, X, R);
  return R;
}

// bar x lou

void barlou(const bar & X, const lou & A, bar & R) {
  if(A.nrow != X.n || A.ncol != R.n) {
    Rcerr << "dim mismatch in barlou\n";
    return;
  }
  R.zeros();
  size_t off = 0;
  for(size_t j = 0; j < A.ncol; j++) 
    for(size_t i = 0; i < A.nrow; i++) 
      R.data[j] += A.data[off++]*X.data[i]; 
}

bar barlou(const bar & X, const lou & A) {
  bar R(A.ncol);
  barlou(X, A, R);
  return R;
}

// bar x bar (produit scalaire)
double barbar(const bar & X, const bar & Y) {
  if(X.n != Y.n) {
    Rcerr << "dim mismatch in barbar\n";
  }
  double r = 0;
  for(size_t i = 0; i < X.n; i++) r += X.data[i]*Y.data[i];
  return r;
}

// lou x lou
void loulou(const lou & A, const lou & B, lou & C) {
  if(A.nrow != C.nrow || A.ncol != B.nrow || B.ncol != C.ncol) {
    Rcerr << "dim mismatch in barlou\n";
    return;
  }
  C.zeros();
  size_t b_off = 0;
  for(size_t j = 0 ; j < C.ncol; j++) {
    size_t a_off = 0;
    for(size_t k = 0 ; k < A.ncol; k++) {
      size_t c_off = j * C.nrow;
      double b_kj = B.data[b_off++];
      for(size_t i = 0 ; i < C.nrow; i++) 
        C.data[c_off++] += A.data[a_off++]*b_kj;
    }
  }
}




/* conversion vers objets R */

NumericMatrix as_r(lou & A) {
  NumericMatrix R(A.nrow, A.ncol);
  double * pR = &R(0,0);
  for(size_t i = 0; i < A.nrow*A.ncol; i++) pR[i] = A.data[i];
  return R;
}

NumericVector as_r(bar & A) {
  NumericVector R(A.n);
  double * pR = &R[0];
  for(size_t i = 0 ; i < A.n; i++) pR[i] = A.data[i];
  return R;
}












