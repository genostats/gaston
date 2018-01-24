#include <Rcpp.h>
#include "Parallel.h"
#include <iostream>
#include <ctime>
#include "matrix4.h"
#include "loubar.h"

using namespace Rcpp;
using namespace Parallel;


struct paraPro_p : public Worker {
  const matrix4 & A;
  const std::vector<double> p;
  const size_t ncol;
  const size_t true_ncol;
  const size_t r;
  const double * v;
  
  //output
  double * Av;

  //constructeur
  paraPro_p(matrix4 & A, std::vector<double> p, size_t r, double * v, double * Av)
           : A(A), p(p), ncol(A.ncol), true_ncol(A.true_ncol), r(r), v(v), Av(Av) {}

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
          for(int ss = 0; ss < 4 && 4*j + ss < ncol; ss++) {
            Av[A.nrow*c+i] += v[k++]*gg[x&3];
            x >>= 2;
          }
        }
      }
    }
  }
};

// produit par des vecteurs de dim p_A->ncol = le nombre d'individus
// c'est pour un produit pour retrouver [par exemple] un loading à partir des PC
// les vecteurs sont supposés en colonne [transposer si besoin avant d'appeler]
// le résultat est en colonne (un loading par colonne)
// [[Rcpp::export]]
NumericMatrix m4_pc_to_loading_p(XPtr<matrix4> p_A, const std::vector<double> & p, NumericMatrix v) {

  int n = p_A->nrow; // nb snps
  int m = p_A->ncol; // nb inds
  if(m != v.nrow()) Rf_error("Dimensions mismatch");
  int r = v.ncol();

  NumericMatrix R(n,r);
  paraPro_p X(*p_A, p, r, v.begin(), R.begin());
  parallelFor(0, n, X, 100);
  return R;
}



RcppExport SEXP gg_m4_pc_to_loading_p(SEXP p_ASEXP, SEXP pSEXP, SEXP vSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< XPtr<matrix4> >::type p_A(p_ASEXP );
        Rcpp::traits::input_parameter< const std::vector<double>& >::type p(pSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type v(vSEXP );
        NumericMatrix __result = m4_pc_to_loading_p(p_A, p, v);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}

