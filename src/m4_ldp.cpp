#include <Rcpp.h>
#include <iostream>
#include <ctime>
#include "matrix4.h"
#include "loubar.h"
#include "m4_ld.h"

using namespace Rcpp;

// la signature les différencie des autres !
void LD_col(matrix4 & A, bar & p, size_t c1, size_t c2, lou & LD) {
  const size_t n = c2-c1+1;
  if(n != LD.nrow || n != LD.ncol) {
    Rcout << "dim mismatch in LD_col (lou)\n" ;
    return;
  }

  for(size_t i1 = 0; i1 < n; i1++) {
    size_t x1 = c1+i1;
    double mu1 = 2*p.data[x1];       
    size_t off = i1*n;
    for(size_t i2 = 0; i2 <= i1; i2++) {
      size_t x2 = c1+i2;
      double mu2 = 2*p.data[x2];       
      double v=2*sqrt(p.data[x1]*(1-p.data[x1])*p.data[x2]*(1-p.data[x2]));
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
void LD_col_0(matrix4 & A, bar & p, int c1, int c2, int d1, int d2, lou & LD) {
  const int n = c2-c1+1;
  const int m = d2-d1+1;
  if(n != LD.nrow || m != LD.ncol) {
    Rcout << "dim mismatch in LD_col_0 (lou)\n" ;
    return;
  }

  for(int i2 = 0; i2 < m; i2++) {
    int x2 = d1+i2;
    double mu2 = 2*p.data[x2];
    size_t off = LD.nrow*i2;  
    for(int i1 = 0; i1 < n; i1++) {
      int x1 = c1+i1;
      double mu1 = 2*p.data[x1];       
      double v = 2*sqrt(p.data[x1]*(1-p.data[x1])*p.data[x2]*(1-p.data[x2]));
      LD.data[off++] = LD_colxx(A, mu1, mu2, v, x1, x2);
    }
  } 
}

// c1 < d1 < c2 < d2
// bcp plus clair avec la boucle sur x1 x2, réécrire les précédentes comme ça ?
void LD_col_1(matrix4 & A, bar & p, int c1, int c2, int d1, int d2, lou & LD) {
  if(c2-c1+1 != LD.nrow || d2-d1+1 != LD.ncol) {
    Rcout << "dim mismatch in LD_col_1 (lou)\n" ;
    return;
  }

  for(int x2 = d1; x2 <= d2; x2++) {
    double mu2 = 2*p.data[x2];
    size_t off = (x2-d1)*LD.nrow;
    for(int x1 = c1; x1 < d1; x1++) {
      double mu1 = 2*p.data[x1];
      double v=2*sqrt(p.data[x1]*(1-p.data[x1])*p.data[x2]*(1-p.data[x2]));
      LD.data[off++] = LD_colxx(A, mu1, mu2, v, x1, x2);
    }
  }
  for(int x2 = d1; x2 <= c2; x2++) {
    double mu2 = 2*p.data[x2];
    size_t off = (x2-d1)*LD.nrow + d1 - c1;
    for(int x1 = d1; x1 <= x2; x1++) {
      double mu1 = 2*p.data[x1];
      double v=2*sqrt(p.data[x1]*(1-p.data[x1])*p.data[x2]*(1-p.data[x2]));
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
    double mu2 = 2*p.data[x2];
    size_t off = (x2-d1)*LD.nrow + d1 - c1;
    for(int x1 = d1; x1 <= c2; x1++) {
      double mu1 = 2*p.data[x1];
      double v=2*sqrt(p.data[x1]*(1-p.data[x1])*p.data[x2]*(1-p.data[x2]));
      LD.data[off++] = LD_colxx(A, mu1, mu2, v, x1, x2);
    }
  }

}


// d1 < c1 < d2 < c2
void LD_col_2(matrix4 & A, bar & p, int c1, int c2, int d1, int d2, lou & LD) {
  if(c2-c1+1 != LD.nrow || d2-d1+1 != LD.ncol) {
    Rcout << "dim mismatch in LD_col_2 (lou)\n" ;
    return;
  }

  for(int x2 = d1; x2 < c1; x2++) {
    double mu2 = 2*p.data[x2];
    size_t off = (x2-d1)*LD.nrow;
    for(int x1 = c1; x1 <= c2; x1++) {
      double mu1 = 2*p.data[x1];
      double v=2*sqrt(p.data[x1]*(1-p.data[x1])*p.data[x2]*(1-p.data[x2]));
      LD.data[off++] = LD_colxx(A, mu1, mu2, v, x1, x2);
    }
  }
  for(int x2 = c1; x2 <= d2; x2++) {
    double mu2 = 2*p.data[x2];
    size_t off = (x2-d1)*LD.nrow;
    for(int x1 = c1; x1 <= x2; x1++) {
      double mu1 = 2*p.data[x1];
      double v=2*sqrt(p.data[x1]*(1-p.data[x1])*p.data[x2]*(1-p.data[x2]));
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
    double mu2 = 2*p.data[x2];
    size_t off = (d2+1-c1) + (x2-d1)*LD.nrow;
    for(int x1 = d2+1; x1 <= c2; x1++) {
      double mu1 = 2*p.data[x1];
      double v=2*sqrt(p.data[x1]*(1-p.data[x1])*p.data[x2]*(1-p.data[x2]));
      LD.data[off++] = LD_colxx(A, mu1, mu2, v, x1, x2);
    }
  }
}

// d1 < c1 < c2 < d2
void LD_col_3(matrix4 & A, bar & p, int c1, int c2, int d1, int d2, lou & LD) {
  if(c2-c1+1 != LD.nrow || d2-d1+1 != LD.ncol) {
    Rcout << "dim mismatch in LD_col_3 (lou)\n" ;
    return;
  }

  for(int x2 = d1; x2 < c1; x2++) {
    double mu2 = 2*p.data[x2];
    size_t off = (x2-d1)*LD.nrow;
    for(int x1 = c1; x1 <= c2; x1++) {
      double mu1 = 2*p.data[x1];
      double v=2*sqrt(p.data[x1]*(1-p.data[x1])*p.data[x2]*(1-p.data[x2]));
      LD.data[off++] = LD_colxx(A, mu1, mu2, v, x1, x2);
    }
  }
  for(int x2 = c1; x2 <= c2; x2++) {
    double mu2 = 2*p.data[x2];
    size_t off = (x2-d1)*LD.nrow;
    for(int x1 = c1; x1 <= x2; x1++) {
      double mu1 = 2*p.data[x1];
      double v=2*sqrt(p.data[x1]*(1-p.data[x1])*p.data[x2]*(1-p.data[x2]));
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
    double mu2 = 2*p.data[x2];
    size_t off = (x2-d1)*LD.nrow;
    for(int x1 = c1; x1 <= c2; x1++) {
      double mu1 = 2*p.data[x1];
      double v=2*sqrt(p.data[x1]*(1-p.data[x1])*p.data[x2]*(1-p.data[x2]));
      LD.data[off++] = LD_colxx(A, mu1, mu2, v, x1, x2);
    }
  }
}


// c1 < d1 < d2 < c2 
void LD_col_4(matrix4 & A, bar & p, int c1, int c2, int d1, int d2, lou & LD) {
  if(c2-c1+1 != LD.nrow || d2-d1+1 != LD.ncol) {
    Rcout << "dim mismatch in LD_col_4 (lou)\n" ;
    return;
  }

  for(int x2 = d1; x2 <= d2; x2++) {
    double mu2 = 2*p.data[x2];
    size_t off = (x2-d1)*LD.nrow;
    for(int x1 = c1; x1 < d1; x1++) {
      double mu1 = 2*p.data[x1];
      double v=2*sqrt(p.data[x1]*(1-p.data[x1])*p.data[x2]*(1-p.data[x2]));
      LD.data[off++] = LD_colxx(A, mu1, mu2, v, x1, x2);
    }
  }
  for(int x2 = d1; x2 <= d2; x2++) {
    double mu2 = 2*p.data[x2];
    size_t off = (d1-c1) + (x2-d1)*LD.nrow;
    for(int x1 = d1; x1 <= x2; x1++) {
      double mu1 = 2*p.data[x1];
      double v=2*sqrt(p.data[x1]*(1-p.data[x1])*p.data[x2]*(1-p.data[x2]));
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
    double mu2 = 2*p.data[x2];
    size_t off = (d2+1-c1) + (x2-d1)*LD.nrow;
    for(int x1 = d2+1; x1 <= c2; x1++) {
      double mu1 = 2*p.data[x1];
      double v=2*sqrt(p.data[x1]*(1-p.data[x1])*p.data[x2]*(1-p.data[x2]));
      LD.data[off++] = LD_colxx(A, mu1, mu2, v, x1, x2);
    }
  }
}


// Choix de la bonne fonction
void LD_chunk(matrix4 & A, bar & p, int c1, int c2, int d1, int d2, lou & LD) {
  if(c2-c1+1 != LD.nrow || d2-d1+1 != LD.ncol) {
    Rcout << "dim mismatch in LD_chunk (lou)\n" ;
    return;
  }

  if(c2 <= d1 || d2 <= c1)
    LD_col_0(A, p, c1, c2, d1, d2, LD);
  else if(c1 <= d1 && c2 <= d2)
    LD_col_1(A, p, c1, c2, d1, d2, LD);
  else if(d1 <= c1 && d2 <= c2)
    LD_col_2(A, p, c1, c2, d1, d2, LD);
  else if(d1 <= c1 && c2 <= d2)
    LD_col_3(A, p, c1, c2, d1, d2, LD);
  else if(c1 <= d1 && d2 <= c2)
    LD_col_4(A, p, c1, c2, d1, d2, LD);

}

// [[Rcpp::export]]
NumericMatrix LD_p(XPtr<matrix4> p_A, NumericVector p, int c1, int c2) {
  bar p_(p, Clone);
  NumericMatrix LD(c2-c1+1, c2-c1+1);
  lou LD_(LD, Clone);
  LD_col(*p_A, p_, c1, c2, LD_);
  return LD;
}

// [[Rcpp::export]]
NumericMatrix LD_chunk_p(XPtr<matrix4> p_A, NumericVector p, int c1, int c2, int d1, int d2) {
  bar p_(p, Clone);
  NumericMatrix LD(c2-c1+1, d2-d1+1);
  lou LD_(LD, Clone);
  LD_chunk(*p_A, p_, c1, c2, d1, d2, LD_);
  return LD;
}

RcppExport SEXP gg_LD_p(SEXP p_ASEXP, SEXP pSEXP, SEXP c1SEXP, SEXP c2SEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< XPtr<matrix4> >::type p_A(p_ASEXP );
        Rcpp::traits::input_parameter< NumericVector >::type p(pSEXP );
        Rcpp::traits::input_parameter< int >::type c1(c1SEXP );
        Rcpp::traits::input_parameter< int >::type c2(c2SEXP );
        NumericMatrix __result = LD_p(p_A, p, c1, c2);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}

RcppExport SEXP gg_LD_chunk_p(SEXP p_ASEXP, SEXP pSEXP, SEXP c1SEXP, SEXP c2SEXP, SEXP d1SEXP, SEXP d2SEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< XPtr<matrix4> >::type p_A(p_ASEXP );
        Rcpp::traits::input_parameter< NumericVector >::type p(pSEXP );
        Rcpp::traits::input_parameter< int >::type c1(c1SEXP );
        Rcpp::traits::input_parameter< int >::type c2(c2SEXP );
        Rcpp::traits::input_parameter< int >::type d1(d1SEXP );
        Rcpp::traits::input_parameter< int >::type d2(d2SEXP );
        NumericMatrix __result = LD_chunk_p(p_A, p, c1, c2, d1, d2);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}

