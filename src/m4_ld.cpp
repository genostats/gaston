#include <Rcpp.h>
#include "Parallel.h"
#include <iostream>
#include <ctime>
#include "matrix4.h"
#include "loubar.h"
#include "m4_ld.h"

using namespace Rcpp;
using namespace Parallel;

/* ************
 * calcul du LD par un crossproduct sur une fenêtre des génotypes
 */

double LD_colxx(matrix4 & A, double mu1, double mu2, double v, size_t x1, size_t x2) {
  double LD = 0;
  double gg[16];
  gg[3] = gg[7] = gg[11] = gg[12] = gg[13] = gg[14] = gg[15] = 0;
  gg[0] = (-mu1)*(-mu2);
  gg[1] = (-mu1)*(1.-mu2);
  gg[2] = (-mu1)*(2.-mu2);

  gg[4] = (1.-mu1)*(-mu2);
  gg[5] = (1.-mu1)*(1.-mu2);
  gg[6] = (1.-mu1)*(2.-mu2);

  gg[8] = (2.-mu1)*(-mu2);
  gg[9] = (2.-mu1)*(1.-mu2);
  gg[10]= (2.-mu1)*(2.-mu2);

  for(size_t i = 0; i < A.true_ncol; i++) {
    uint8_t g1 = A.data[x1][i];
    uint8_t g2 = A.data[x2][i];
    for(int ss = 0; ss < 4; ss++) {
      LD += gg[ ((g1&3)*4) + (g2&3) ];
      g1 >>= 2;
      g2 >>= 2;
    }
  }
  return LD/(v*(A.ncol-1));   
}


// la version parallèle est plus lente (en tout cas sur mon portable)... (c'est le join ?)
struct paraLD : public Worker {
  // input
  const uint8_t * d1;
  const uint8_t * d2;
  double * gg;

  // output
  double LD;

  paraLD(matrix4 & A, double * gg, int x1, int x2) : d1(A.data[x1]), d2(A.data[x2]), gg(gg), LD(0) {}
  paraLD(paraLD & Q, Split) : d1(Q.d1), d2(Q.d2), gg(Q.gg), LD(0) {}

  void operator()(size_t beg, size_t end) {
    for(size_t i = beg; i < end; i++) {
      uint8_t g1 = d1[i];
      uint8_t g2 = d2[i];
      for(int ss = 0; ss < 4; ss++) {
        LD += gg[ ((g1&3)*4) + (g2&3) ];
        g1 >>= 2;
        g2 >>= 2;
      }
    }   
  }
  void join(const paraLD & Q) { LD += Q.LD; }
};


inline double LD_colxx1(matrix4 & A, double mu1, double mu2, double v, size_t x1, size_t x2) {
  double gg[16];
  gg[3] = gg[7] = gg[11] = gg[12] = gg[13] = gg[14] = gg[15] = 0;
  gg[0] = (-mu1)*(-mu2);
  gg[1] = (-mu1)*(1.-mu2);
  gg[2] = (-mu1)*(2.-mu2);

  gg[4] = (1.-mu1)*(-mu2);
  gg[5] = (1.-mu1)*(1.-mu2);
  gg[6] = (1.-mu1)*(2.-mu2);

  gg[8] = (2.-mu1)*(-mu2);
  gg[9] = (2.-mu1)*(1.-mu2);
  gg[10]= (2.-mu1)*(2.-mu2);

  paraLD X(A, gg, x1, x2);
  parallelReduce(0,A.true_ncol,X,1);
  return X.LD/(v*(A.ncol-1));  
}








/**************************************************************
 *
 *   on refait tout avec une signature lou / bar
 *
 **************************************************************/

void LD_col(matrix4 & A, bar & mu, bar & sd, size_t c1, size_t c2, lou & LD) {
  const size_t n = c2-c1+1;
  if(n != LD.nrow || n != LD.ncol) {
    Rcout << "dim mismatch in LD_col (lou)\n" ;
    return;
  }

  for(size_t i1 = 0; i1 < n; i1++) {
    size_t x1 = c1+i1;
    double mu1 = mu.data[x1];       
    size_t off = i1*n;
    for(size_t i2 = 0; i2 <= i1; i2++) {
      size_t x2 = c1+i2;
      double mu2 = mu.data[x2];       
      double v=sd.data[x1]*sd.data[x2];
      LD.data[off++] = LD_colxx(A, mu1, mu2, v, x1, x2);
    }
  } 

  // symmetriser
  for(size_t i1 = 0; i1 < n; i1++) {
    for(size_t i2 = 0; i2 < i1; i2++) {
       LD.data[i1 + n*i2] = LD.data[i2 + n*i1];
    }
  }
}


// Intervalles c1 c2 et d1 d2 disjoints
void LD_col_0(matrix4 & A, bar & mu, bar & sd, int c1, int c2, int d1, int d2, lou & LD) {
  const int n = c2-c1+1;
  const int m = d2-d1+1;
  if(n != LD.nrow || m != LD.ncol) {
    Rcout << "dim mismatch in LD_col_0 (lou)\n" ;
    return;
  }

  for(int i2 = 0; i2 < m; i2++) {
    int x2 = d1+i2;
    double mu2 = mu.data[x2];
    size_t off = LD.nrow*i2;  
    for(int i1 = 0; i1 < n; i1++) {
      int x1 = c1+i1;
      double mu1 = mu.data[x1];       
      double v = sd.data[x1]*sd.data[x2];
      LD.data[off++] = LD_colxx(A, mu1, mu2, v, x1, x2);
    }
  } 
}

// c1 < d1 < c2 < d2
// bcp plus clair avec la boucle sur x1 x2, réécrire les précédentes comme ça ?
void LD_col_1(matrix4 & A, bar & mu, bar & sd, int c1, int c2, int d1, int d2, lou & LD) {
  if(c2-c1+1 != LD.nrow || d2-d1+1 != LD.ncol) {
    Rcout << "dim mismatch in LD_col_1 (lou)\n" ;
    return;
  }

  for(int x2 = d1; x2 <= d2; x2++) {
    double mu2 = mu.data[x2];
    size_t off = (x2-d1)*LD.nrow;
    for(int x1 = c1; x1 < d1; x1++) {
      double mu1 = mu.data[x1];
      double v=sd.data[x1]*sd.data[x2];
      LD.data[off++] = LD_colxx(A, mu1, mu2, v, x1, x2);
    }
  }
  for(int x2 = d1; x2 <= c2; x2++) {
    double mu2 = mu.data[x2];
    size_t off = (x2-d1)*LD.nrow + d1 - c1;
    for(int x1 = d1; x1 <= x2; x1++) {
      double mu1 = mu.data[x1];
      double v=sd.data[x1]*sd.data[x2];
      LD.data[off++] = LD_colxx(A, mu1, mu2, v, x1, x2);
    }
  }
  // symmetriser ce morceau
  for(int x1 = d1; x1 <= c2; x1++) {
    for(int x2 = d1; x2 < x1; x2++) {
      LD.data[ (x1-c1) + LD.nrow*(x2-d1) ] = LD.data[ (x2-c1) + LD.nrow*(x1-d1) ];
    }
  }
  for(int x2 = c2+1; x2 <= d2; x2++) {
    double mu2 = mu.data[x2];
    size_t off = (x2-d1)*LD.nrow + d1 - c1;
    for(int x1 = d1; x1 <= c2; x1++) {
      double mu1 = mu.data[x1];
      double v=sd.data[x1]*sd.data[x2];
      LD.data[off++] = LD_colxx(A, mu1, mu2, v, x1, x2);
    }
  }

}


// d1 < c1 < d2 < c2
void LD_col_2(matrix4 & A, bar & mu, bar & sd, int c1, int c2, int d1, int d2, lou & LD) {
  if(c2-c1+1 != LD.nrow || d2-d1+1 != LD.ncol) {
    Rcout << "dim mismatch in LD_col_2 (lou)\n" ;
    return;
  }

  for(int x2 = d1; x2 < c1; x2++) {
    double mu2 = mu.data[x2];
    size_t off = (x2-d1)*LD.nrow;
    for(int x1 = c1; x1 <= c2; x1++) {
      double mu1 = mu.data[x1];
      double v=sd.data[x1]*sd.data[x2];
      LD.data[off++] = LD_colxx(A, mu1, mu2, v, x1, x2);
    }
  }
  for(int x2 = c1; x2 <= d2; x2++) {
    double mu2 = mu.data[x2];
    size_t off = (x2-d1)*LD.nrow;
    for(int x1 = c1; x1 <= x2; x1++) {
      double mu1 = mu.data[x1];
      double v=sd.data[x1]*sd.data[x2];
      LD.data[off++] = LD_colxx(A, mu1, mu2, v, x1, x2);
    }
  }
  // symmetriser ce morceau
  for(int x1 = c1; x1 <= d2; x1++) {
    for(int x2 = c1; x2 < x1; x2++) {
      LD.data[ (x1-c1) + (x2-d1)*LD.nrow ] = LD.data[ (x2-c1) + (x1-d1)*LD.nrow ];
    }
  }
  for(int x2 = c1; x2 <= d2; x2++) {
    double mu2 = mu.data[x2];
    size_t off = (d2+1-c1) + (x2-d1)*LD.nrow;
    for(int x1 = d2+1; x1 <= c2; x1++) {
      double mu1 = mu.data[x1];
      double v=sd.data[x1]*sd.data[x2];
      LD.data[off++] = LD_colxx(A, mu1, mu2, v, x1, x2);
    }
  }
}

// d1 < c1 < c2 < d2
void LD_col_3(matrix4 & A, bar & mu, bar & sd, int c1, int c2, int d1, int d2, lou & LD) {
  if(c2-c1+1 != LD.nrow || d2-d1+1 != LD.ncol) {
    Rcout << "dim mismatch in LD_col_3 (lou)\n" ;
    return;
  }

  for(int x2 = d1; x2 < c1; x2++) {
    double mu2 = mu.data[x2];
    size_t off = (x2-d1)*LD.nrow;
    for(int x1 = c1; x1 <= c2; x1++) {
      double mu1 = mu.data[x1];
      double v=sd.data[x1]*sd.data[x2];
      LD.data[off++] = LD_colxx(A, mu1, mu2, v, x1, x2);
    }
  }
  for(int x2 = c1; x2 <= c2; x2++) {
    double mu2 = mu.data[x2];
    size_t off = (x2-d1)*LD.nrow;
    for(int x1 = c1; x1 <= x2; x1++) {
      double mu1 = mu.data[x1];
      double v=sd.data[x1]*sd.data[x2];
      LD.data[off++] = LD_colxx(A, mu1, mu2, v, x1, x2);
    }
  }
  // symmetriser ce morceau
  for(int x1 = c1; x1 <= c2; x1++) {
    for(int x2 = c1; x2 < x1; x2++) {
      LD.data[ (x1-c1) + (x2-d1)*LD.nrow ] = LD.data[ (x2-c1) + (x1-d1)*LD.nrow ];
    }
  }
  for(int x2 = c2+1; x2 <= d2; x2++) {
    double mu2 = mu.data[x2];
    size_t off = (x2-d1)*LD.nrow;
    for(int x1 = c1; x1 <= c2; x1++) {
      double mu1 = mu.data[x1];
      double v=sd.data[x1]*sd.data[x2];
      LD.data[off++] = LD_colxx(A, mu1, mu2, v, x1, x2);
    }
  }
}


// c1 < d1 < d2 < c2 
void LD_col_4(matrix4 & A, bar & mu, bar & sd, int c1, int c2, int d1, int d2, lou & LD) {
  if(c2-c1+1 != LD.nrow || d2-d1+1 != LD.ncol) {
    Rcout << "dim mismatch in LD_col_4 (lou)\n" ;
    return;
  }

  for(int x2 = d1; x2 <= d2; x2++) {
    double mu2 = mu.data[x2];
    size_t off = (x2-d1)*LD.nrow;
    for(int x1 = c1; x1 < d1; x1++) {
      double mu1 = mu.data[x1];
      double v=sd.data[x1]*sd.data[x2];
      LD.data[off++] = LD_colxx(A, mu1, mu2, v, x1, x2);
    }
  }
  for(int x2 = d1; x2 <= d2; x2++) {
    double mu2 = mu.data[x2];
    size_t off = (d1-c1) + (x2-d1)*LD.nrow;
    for(int x1 = d1; x1 <= x2; x1++) {
      double mu1 = mu.data[x1];
      double v=sd.data[x1]*sd.data[x2];
      LD.data[off++] = LD_colxx(A, mu1, mu2, v, x1, x2);
    }
  }
  // symmetriser ce morceau
  for(int x1 = d1; x1 <= d2; x1++) {
    for(int x2 = d1; x2 < x1; x2++) {
      LD.data[ (x1-c1) + (x2-d1)*LD.nrow ] = LD.data[ (x2-c1) + (x1-d1)*LD.nrow ];
    }
  }
  for(int x2 = d1; x2 <= d2; x2++) {
    double mu2 = mu.data[x2];
    size_t off = (d2+1-c1) + (x2-d1)*LD.nrow;
    for(int x1 = d2+1; x1 <= c2; x1++) {
      double mu1 = mu.data[x1];
      double v=sd.data[x1]*sd.data[x2];
      LD.data[off++] = LD_colxx(A, mu1, mu2, v, x1, x2);
    }
  }
}


// Choix de la bonne fonction
void LD_chunk(matrix4 & A, bar & mu, bar & sd, int c1, int c2, int d1, int d2, lou & LD) {
  if(c2-c1+1 != LD.nrow || d2-d1+1 != LD.ncol) {
    Rcout << "dim mismatch in LD_chunk (lou)\n" ;
    return;
  }

  if(c2 <= d1 || d2 <= c1)
    LD_col_0(A, mu, sd, c1, c2, d1, d2, LD);
  else if(c1 <= d1 && c2 <= d2)
    LD_col_1(A, mu, sd, c1, c2, d1, d2, LD);
  else if(d1 <= c1 && d2 <= c2)
    LD_col_2(A, mu, sd, c1, c2, d1, d2, LD);
  else if(d1 <= c1 && c2 <= d2)
    LD_col_3(A, mu, sd, c1, c2, d1, d2, LD);
  else if(c1 <= d1 && d2 <= c2)
    LD_col_4(A, mu, sd, c1, c2, d1, d2, LD);

}

// [[Rcpp::export]]
NumericMatrix LD(XPtr<matrix4> p_A, NumericVector mu, NumericVector sd, int c1, int c2) {
  bar mu_(mu, Clone);
  bar sd_(sd, Clone);
  NumericMatrix LD(c2-c1+1, c2-c1+1);
  lou LD_(LD, Clone);
  LD_col(*p_A, mu_, sd_, c1, c2, LD_);
  return LD;
}

// [[Rcpp::export]]
NumericMatrix LD_chunk(XPtr<matrix4> p_A, NumericVector mu, NumericVector sd, int c1, int c2, int d1, int d2) {
  bar mu_(mu, Clone);
  bar sd_(sd, Clone);
  NumericMatrix LD(c2-c1+1, d2-d1+1);
  lou LD_(LD, Clone);
  LD_chunk(*p_A, mu_, sd_, c1, c2, d1, d2, LD_);
  return LD;
}



RcppExport SEXP gg_LD(SEXP p_ASEXP, SEXP muSEXP, SEXP sdSEXP, SEXP c1SEXP, SEXP c2SEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< XPtr<matrix4> >::type p_A(p_ASEXP );
        Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type sd(sdSEXP );
        Rcpp::traits::input_parameter< int >::type c1(c1SEXP );
        Rcpp::traits::input_parameter< int >::type c2(c2SEXP );
        NumericMatrix __result = LD(p_A, mu, sd, c1, c2);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}

RcppExport SEXP gg_LD_chunk(SEXP p_ASEXP, SEXP muSEXP, SEXP sdSEXP, SEXP c1SEXP, SEXP c2SEXP, SEXP d1SEXP, SEXP d2SEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< XPtr<matrix4> >::type p_A(p_ASEXP );
        Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type sd(sdSEXP );
        Rcpp::traits::input_parameter< int >::type c1(c1SEXP );
        Rcpp::traits::input_parameter< int >::type c2(c2SEXP );
        Rcpp::traits::input_parameter< int >::type d1(d1SEXP );
        Rcpp::traits::input_parameter< int >::type d2(d2SEXP );
        NumericMatrix __result = LD_chunk(p_A, mu, sd, c1, c2, d1, d2);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}

