#include <Rcpp.h>
#include <iostream>
#include <ctime>
#include "matrix4.h"
#include "loubar.h"

using namespace Rcpp;

// [[Rcpp::export]]
XPtr<matrix4> new_matrix4(int nrow, int ncol) {
  XPtr<matrix4> p_A(new matrix4(nrow, ncol));
  return p_A;
}

// [[Rcpp::export]]
XPtr<matrix4> as_matrix4(NumericMatrix A) {
  XPtr<matrix4> p_x(new matrix4(A));
  return p_x;
}

// [[Rcpp::export]]
XPtr<matrix4> raw_as_matrix4(RawMatrix A) {
  XPtr<matrix4> p_x(new matrix4(A));
  return p_x;
}

// une conversion qui transpose pour un rendu plus naturel 
// et des cheveux plus brillants
// [[Rcpp::export]]
IntegerMatrix m4_as012(XPtr<matrix4> pA) {
  IntegerMatrix X(pA->ncol, pA->nrow);
  for(int j = 0; j < X.ncol(); j++) {
    for(int i = 0; i < pA->true_ncol-1; i++) {
      uint8_t x = pA->data[j][i];
      for(int ss = 0; ss < 4; ss++) {
        X(4*i+ss,j) = ((x&3) != 3)?(x&3):NA_INTEGER;
        x >>= 2;
      }
    }
    int i = pA->true_ncol-1;
    uint8_t x = pA->data[j][i];
    for(int ss = 0; ss < 4 && 4*i+ss < pA->ncol; ss++) {
      X(4*i+ss,j) = ((x&3) != 3)?(x&3):NA_INTEGER;
      x >>= 2;
    }
  }
  return X;
}

// [[Rcpp::export]]
NumericMatrix m4_as_scaled_matrix_p(XPtr<matrix4> pA, NumericVector p) {
  if(p.length() != pA->nrow) 
    Rf_error("Dimension mismatch");
  NumericMatrix X(pA->ncol, pA->nrow);
  for(int j = 0; j < X.ncol(); j++) {
    double gg[4] = {   -2*p[j] /sqrt(2*p[j]*(1-p[j])), 
                     (1-2*p[j])/sqrt(2*p[j]*(1-p[j])), 
                     (2-2*p[j])/sqrt(2*p[j]*(1-p[j])), NA_REAL };
    for(int i = 0; i < pA->true_ncol-1; i++) {
      uint8_t x = pA->data[j][i];
      for(int ss = 0; ss < 4; ss++) {
        X(4*i+ss,j) = gg[x&3]; // ((x&3) != 3)?gg[x&3]:NA_REAL;
        x >>= 2;
      }
    }
    int i = pA->true_ncol-1;
    uint8_t x = pA->data[j][i];
    for(int ss = 0; ss < 4 && 4*i+ss < pA->ncol; ss++) {
      X(4*i+ss,j) = gg[x&3]; // ((x&3) != 3)?gg[x&3]:NA_REAL;
      x >>= 2;
    }
  }
  return X;
}


// [[Rcpp::export]]
NumericMatrix m4_as_scaled_matrix_mu_sigma(XPtr<matrix4> pA, NumericVector mu, NumericVector sigma) {
  if(mu.length() != pA->nrow || sigma.length() != pA->nrow) 
    Rf_error("Dimension mismatch");
  NumericMatrix X(pA->ncol, pA->nrow);
  for(int j = 0; j < X.ncol(); j++) {
    double gg[4] = { -mu[j]/sigma[j], (1-mu[j])/sigma[j], (2-mu[j])/sigma[j], NA_REAL };
    for(int i = 0; i < pA->true_ncol-1; i++) {
      uint8_t x = pA->data[j][i];
      for(int ss = 0; ss < 4; ss++) {
        X(4*i+ss,j) = gg[x&3]; // ((x&3) != 3)?gg[x&3]:NA_REAL;
        x >>= 2;
      }
    }
    int i = pA->true_ncol-1;
    uint8_t x = pA->data[j][i];
    for(int ss = 0; ss < 4 && 4*i+ss < pA->ncol; ss++) {
      X(4*i+ss,j) = gg[x&3]; // ((x&3) != 3)?gg[x&3]:NA_REAL;
      x >>= 2;
    }
  }
  return X;
}



