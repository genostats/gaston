#include <Rcpp.h>
#include "Parallel.h"
#include <iostream>
#include <ctime>
#include "matrix4.h"
#include "loubar.h"

using namespace Rcpp;
using namespace Parallel;


struct paraPro2_p : public Worker {
  const matrix4 & A;
  const std::vector<double> p;
  const size_t ncol;
  const size_t true_ncol;
  const size_t r;
  double * v;
  //output
  double * vA;

 //constructeur
  paraPro2_p(matrix4 & A, std::vector<double> p, size_t r, double * v) 
            : A(A), p(p), ncol(A.ncol), true_ncol(A.true_ncol), r(r), v(v) {
    vA = new double[r*ncol];
    std::fill(vA, vA+r*ncol, 0); 
  }

  //constructeur pour le split
  paraPro2_p(paraPro2_p & Q, Split) : A(Q.A), p(Q.p), ncol(Q.ncol),
           true_ncol(Q.true_ncol), r(Q.r), v(Q.v) {
    vA = new double[r*ncol];
    std::fill(vA, vA+r*ncol, 0);
  }

  // destructeur
  ~paraPro2_p() {
    delete [] vA;
  }

  //worker
  void operator()(size_t beg, size_t end) {
    double gg[4];
    gg[3] = 0;
    for(size_t i = beg; i < end; i++) {
      double mu_ = 2*p[i];
      double sd_ = (mu_ == 0 || mu_ == 2)?1:sqrt(mu_*(1-mu_/2));
      gg[0] = -mu_/sd_;
      gg[1] = (1-mu_)/sd_;
      gg[2] = (2-mu_)/sd_;
      for(size_t c = 0; c < r; c++) {
        size_t k = c*ncol;
        for(size_t j = 0; j < true_ncol; j++) {
          uint8_t x = A.data[i][j];
          for(int ss = 0; ss < 4 && (4*j + ss) < ncol; ss++) {
            vA[k++] += v[A.nrow*c+i]*gg[x&3];
            x >>= 2;
          }
        }
      }
    }
  }

  // join
  void join(const paraPro2_p & Q) {
    std::transform(vA, vA + r*ncol, Q.vA, vA, std::plus<double>());
    // autrement dit : vA += vA.K;
  }

};


//[[Rcpp::export]]
NumericMatrix m4_loading_to_pc_p(XPtr<matrix4> p_A, const std::vector<double> & p, NumericMatrix & v) {
  int n = p_A->nrow; // nb snps
  int m = p_A->ncol; // nb inds
  if(n != v.nrow()) stop("Dimensions mismatch");
  int r = v.ncol();

  paraPro2_p X(*p_A, p, r, v.begin());
  parallelReduce(0, n, X);

  NumericMatrix R(m,r);
  std::copy(X.vA, X.vA+m*r, R.begin());

  return R;
}

RcppExport SEXP gg_m4_loading_to_pc_p(SEXP p_ASEXP, SEXP pSEXP, SEXP vSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< XPtr<matrix4> >::type p_A(p_ASEXP );
        Rcpp::traits::input_parameter< const std::vector<double>& >::type p(pSEXP );
        Rcpp::traits::input_parameter< NumericMatrix& >::type v(vSEXP );
        NumericMatrix __result = m4_loading_to_pc_p(p_A, p, v);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}

