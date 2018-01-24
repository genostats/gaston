#include <Rcpp.h>
#include <iostream>
#include <ctime>
#include "matrix4.h"
#include "loubar.h"

using namespace Rcpp;

/*****************************************************************/

matrix4 transposexx(matrix4 & A) {
  matrix4 B(A.ncol, A.nrow);
  for(size_t j = 0; j < A.nrow; j++) {
    for(size_t i = 0; i < A.true_ncol-1; i++) {
      uint8_t x = A.data[j][i];
      for(int ss = 0; ss < 4; ss++) {
        B(4*i+ss,j) = (x&3);
        x >>= 2;
      }
    }
    size_t i = A.true_ncol-1;
    uint8_t x = A.data[j][i];
    for(int ss = 0; 4*i + ss < A.ncol; ss++) {
      B(4*i+ss,j) = (x&3);
      x >>= 2;
    }
  }
  return B;
}

// [[Rcpp::export]]
XPtr<matrix4> transpose_m4(XPtr<matrix4> p_A) {
  XPtr<matrix4> p_x(new matrix4(transposexx(*p_A)));
  // *p_x = transposexx(*p_A);
  return p_x;
}
