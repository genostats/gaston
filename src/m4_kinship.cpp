#include <Rcpp.h>
#include "Parallel.h"
#include <iostream>
#include "matrix4.h"
#include "m4_kinship_type.h"

using namespace Rcpp;
using namespace Parallel;


// ************* [on ne symmétrise pas] ***********

struct paraKin : public Worker {
  // input and others
  uint8_t ** data;
  const size_t ncol;
  const size_t true_ncol;
  const double * mu;
  const double * w;
  const size_t sizeK;

  // output
  Ktype * K;
  
  // constructeurs
  paraKin(uint8_t ** data, const size_t ncol, const size_t true_ncol, const double * mu, const double * w) : 
        data(data), ncol(ncol), true_ncol(true_ncol), mu(mu), w(w), sizeK((4*true_ncol)*(4*true_ncol+1)/2) { 
          K = new Ktype[sizeK];  // K is padded to a multiple of 4...
          std::fill(K, K+sizeK, 0);
        }
  paraKin(paraKin & Q, Split) : data(Q.data), ncol(Q.ncol), true_ncol(Q.true_ncol), mu(Q.mu), w(Q.w), sizeK(Q.sizeK) {
          K = new Ktype[sizeK];  // K is padded to a multiple of 4...
          std::fill(K, K+sizeK, 0);
        }

  // destructeur
  ~paraKin() { 
          delete [] K; 
  }

  // worker !
  void operator()(size_t beg, size_t end) {
    Ktype H[4];
    Ktype h0[32];
    Ktype h1[32];
    H[3] = 0;

    for(size_t i = beg; i < end; i++) {
      Ktype w_ = (Ktype) w[i]; 
      if(w_ == 0) continue;
      Ktype mu_ = (Ktype) mu[i];
      Ktype v0 = -mu_*w_;
      Ktype v1 = (1-mu_)*w_;
      Ktype v2 = (2-mu_)*w_;
      H[0] = v0;
      H[1] = v1;
      H[2] = v2;
      
      size_t k = 0;
      for(size_t a1 = 0; a1 < 4; a1++) {
        for(size_t a2 = 0; a2 < 4; a2++) {
          h0[k++] = H[a2];
          h0[k++] = H[a1];
        }
      }

      const uint8_t * dd = data[i];
      
      k = 0;
      for(size_t j1 = 0; j1 < true_ncol; j1++) {
        uint8_t x1 = dd[j1];
        for(unsigned int ss1 = 0; (ss1 < 4); ss1++) {
          for(size_t a = 0; a < 32; a++) h1[a] = H[x1&3]*h0[a];
          for(size_t j2 = 0; j2 < j1; j2++) {
            uint8_t x2 = dd[j2];
            K[k++] += h1[ (x2&15)<<1 ];
            K[k++] += h1[ ((x2&15)<<1)+1 ];
            K[k++] += h1[ (x2>>4)<<1 ];
            K[k++] += h1[ ((x2>>4)<<1)+1 ];
          }
          size_t j2 = j1;
          uint8_t x2 = dd[j2];
          for(unsigned int ss2 = 0; ss2 <= ss1; ss2++) {
            K[k++] += H[x1&3]*H[x2&3]; // h1[(x2&3)<<1];
            x2 >>= 2;
          } 
          x1 >>= 2; 
        }
      } 
    }
  }

  // recoller
  void join(const paraKin & Q) {
    std::transform(K, K + sizeK, Q.K, K, std::plus<Ktype>());
    // autrement dit : K += Q.K;
  }

};


NumericMatrix Kinship(XPtr<matrix4> p_A, const std::vector<double> & mu, const std::vector<double> & w, int chunk) {

  if(mu.size() != p_A->nrow || w.size() != p_A->nrow) 
    stop("Dimensions mismatch");

  paraKin X(p_A->data, p_A->ncol, p_A->true_ncol, &mu[0], &w[0]);
  parallelReduce(0, p_A->nrow, X, chunk);

  NumericMatrix Y(p_A->ncol,p_A->ncol);
  size_t k = 0;
  for(size_t i = 0; i < p_A->ncol; i++) {
    for(size_t j = 0; j <= i; j++) {
      Y(j,i) = (double) X.K[k++];
    }
  }

  // symmetriser
  k = 0;
  for(size_t i = 0; i < p_A->ncol; i++) {
    for(size_t j = 0; j <= i; j++) {
      Y(i,j) = X.K[k++]; // ou Y(j,i);
    }
  }
  
  return Y;
}

// [[Rcpp::export]]
NumericMatrix Kinship_w(XPtr<matrix4> p_A, const std::vector<double> & mu, const std::vector<double> & w, LogicalVector snps, int chunk) {
  int nb_snps = sum(snps);

  if(snps.length() != p_A->nrow || mu.size() != nb_snps || w.size() != nb_snps) 
    stop("Dimensions mismatch");

  uint8_t ** data = new uint8_t * [nb_snps];
  size_t k = 0;
  for(size_t i = 0; i < p_A->nrow; i++) {
    if(snps[i]) data[k++] = p_A->data[i];
  }

  paraKin X(data, p_A->ncol, p_A->true_ncol, &mu[0], &w[0]);
  parallelReduce(0, nb_snps, X, chunk);

  delete [] data;

  NumericMatrix Y(p_A->ncol,p_A->ncol);
  k = 0;
  for(size_t i = 0; i < p_A->ncol; i++) {
    for(size_t j = 0; j <= i; j++) {
      Y(j,i) = (double) X.K[k++];
    }
  }

  // symmetriser
  k = 0;
  for(size_t i = 0; i < p_A->ncol; i++) {
    for(size_t j = 0; j <= i; j++) {
      Y(i,j) = X.K[k++]; // ou Y(j,i);
    }
  }
  
  return Y;
}

RcppExport SEXP gg_Kinship(SEXP p_ASEXP, SEXP muSEXP, SEXP wSEXP, SEXP chunkSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< XPtr<matrix4> >::type p_A(p_ASEXP );
        Rcpp::traits::input_parameter< const std::vector<double>& >::type mu(muSEXP );
        Rcpp::traits::input_parameter< const std::vector<double>& >::type w(wSEXP );
        Rcpp::traits::input_parameter< int >::type chunk(chunkSEXP );
        NumericMatrix __result = Kinship(p_A, mu, w, chunk);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}

RcppExport SEXP gg_Kinship_w(SEXP p_ASEXP, SEXP muSEXP, SEXP wSEXP, SEXP snpsSEXP, SEXP chunkSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< XPtr<matrix4> >::type p_A(p_ASEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type w(wSEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type snps(snpsSEXP);
    Rcpp::traits::input_parameter< int >::type chunk(chunkSEXP);
    __result = Rcpp::wrap(Kinship_w(p_A, mu, w, snps, chunk));
    return __result;
END_RCPP
}

