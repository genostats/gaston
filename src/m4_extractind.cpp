#include <Rcpp.h>
#include <iostream>
#include <ctime>
#include "matrix4.h"

using namespace Rcpp;

// [[Rcpp::export]]
XPtr<matrix4> extract_inds_bool(XPtr<matrix4> pA, LogicalVector w) {
  size_t ncol = sum(w);
  if(w.length() != pA->ncol) 
    Rf_error("Length of logical vector doesn't match number of individuals");

  XPtr<matrix4> pB(new matrix4(pA->nrow, ncol));
  for(size_t i = 0; i < pA->nrow; i++){
    size_t k = 0;
    for(size_t j = 0; j < pA->ncol; j++) {
      if(w(j)) {
        (*pB)(i, k) = (*pA)(i,j);
        k++;
      }
    }
  }
  return pB;
}

// [[Rcpp::export]]
XPtr<matrix4> extract_inds_indices(XPtr<matrix4> pA, IntegerVector w) {
  size_t ncol = w.length();
  XPtr<matrix4> pB(new matrix4(pA->nrow, ncol));
  if(is_true(any(w > pA-> ncol)))
     stop("Index out of range"); 
  for(size_t i =0; i < pA->nrow; i++){
    for(size_t j = 0; j < ncol; j++) {
      if(w(j) < 1) 
        (*pB)(i, j) = 3;
      else
        (*pB)(i, j) = (*pA)(i,w(j)-1);
    }
  }
  return pB;
}
