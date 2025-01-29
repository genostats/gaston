#include <Rcpp.h>
#include <iostream>
#include <ctime>
#include "matrix4.h"
#include "loubar.h"

using namespace Rcpp;

/***********************************************************/

// [[Rcpp::export]]
int nsnps(XPtr<matrix4> p_A) {
  return p_A->nrow;
}

// [[Rcpp::export]]
int ninds(XPtr<matrix4> p_A) {
  return p_A->ncol;
}


