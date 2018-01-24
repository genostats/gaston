#include <Rcpp.h>
#include <iostream>
#include <ctime>
#include "matrix4.h"
#include "loubar.h"

using namespace Rcpp;

/***********************************************************/

// [[Rcpp::export]]
void cat_matrix4(XPtr<matrix4> p_A) {
  Rcout << *p_A ;
}

// [[Rcpp::export]]
void fill_line(XPtr<matrix4> p_A, size_t li, NumericVector w) {
  p_A->fill_line(li, w);
}

