#include <Rcpp.h>
#include <iostream>
#include <ctime>
#include "matrix4.h"
#include "loubar.h"

using namespace Rcpp;


// [[Rcpp::export]]
XPtr<matrix4> extract_snps_bool(XPtr<matrix4> pA, LogicalVector w) {
  size_t nrow = sum(w);
  if(w.length() != pA->nrow) 
    Rf_error("Length of logical vector doesn't match number of SNPs");

  XPtr<matrix4> pB(new matrix4(nrow, pA->ncol));
  size_t k = 0;
  for(size_t i =0; i < pA->nrow; i++){
    if(w(i)) {
      std::copy(pA->data[i], pA->data[i]+pA->true_ncol, pB->data[k]);
      k++;
    }
  }
  return pB;
}

// [[Rcpp::export]]
XPtr<matrix4> extract_snps_indices(XPtr<matrix4> pA, IntegerVector w) {
  size_t nrow = w.length();
  XPtr<matrix4> pB(new matrix4(nrow, pA->ncol));
  for(size_t i =0; i < nrow; i++){
    if(w(i) < 1 || w(i) > pA->nrow)
      Rf_error("Index out of range");
    std::copy(pA->data[w(i)-1], pA->data[w(i)-1]+pA->true_ncol, pB->data[i]);
  }
  return pB;
}


