#include <Rcpp.h>
#include "Parallel.h"
#include <iostream>
#include "matrix4.h"
#include "m4_kinship_type.h"

using namespace Rcpp;
using namespace Parallel;


// ************* [on ne symmétrise pas] ***********

struct paraKin_p : public Worker {
  // input and others
  uint8_t ** data;
  const size_t ncol;
  const size_t true_ncol;
  const Ktype n;                  // should be (nb_snps - 1)
  const double * p;
  const size_t sizeK;
  const bool dominance;

  // output
  Ktype * K;
  
  // constructeurs
  paraKin_p(uint8_t ** data, const size_t ncol, const size_t true_ncol, const Ktype n, const double * p, const bool dominance) :
          data(data), ncol(ncol), true_ncol(true_ncol), n(n), p(p), sizeK((4*true_ncol)*(4*true_ncol+1)/2) , dominance(dominance) { 
          K = new Ktype[sizeK];  // K is padded to a multiple of 4...
          std::fill(K, K+sizeK, 0);
        }
  paraKin_p(paraKin_p & Q, Split) : data(Q.data), ncol(Q.ncol), true_ncol(Q.true_ncol), n(Q.n), p(Q.p), sizeK(Q.sizeK), dominance(Q.dominance) {
          K = new Ktype[sizeK];  
          std::fill(K, K+sizeK, 0);
        }

  // destructeur
  ~paraKin_p() { 
          delete [] K; 
  }

  // worker !
  void operator()(size_t beg, size_t end) {
    Ktype H[4];
    Ktype h0[32];
    Ktype h1[32];
    H[3] = 0;

    for(size_t i = beg; i < end; i++) {
      Ktype p_ = (Ktype) p[i];   // freq de l'allele 1 (noté q le plus souvent !!)
      if(p_ == 0 || p_ == 1) continue;
      Ktype v0, v1, v2;
      if(dominance) {
        v0 = p_/(1-p_);
        v1 = -1;
        v2 = (1-p_)/p_; 
      } else {
        Ktype w_ = 1/sqrt(2*p_*(1-p_)); 
        v0 = -2*p_*w_;
        v1 = (1-2*p_)*w_;
        v2 = (2-2*p_)*w_;
      }
      // on divise par le nombre de SNPs
      H[0] = v0/sqrt(n);
      H[1] = v1/sqrt(n);
      H[2] = v2/sqrt(n);
      
      size_t k = 0;
      for(size_t a1 = 0; a1 < 4; a1++) {
        for(size_t a2 = 0; a2 < 4; a2++) {
          h0[k++] = H[a2];
          h0[k++] = H[a1];
        }
      }

      uint8_t * dd = data[i];
      
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
  void join(const paraKin_p & Q) {
    std::transform(K, K + sizeK, Q.K, K, std::plus<Ktype>());
    // autrement dit : K += Q.K;
  }

};

//[[Rcpp::export]]
NumericMatrix Kinship_pw(XPtr<matrix4> p_A, const std::vector<double> & p, LogicalVector snps, bool dominance, int chunk) {
  int nb_snps = sum(snps);

  if(snps.length() != p_A->nrow || p.size() != nb_snps)
    stop("Dimensions mismatch");

  uint8_t ** data = new uint8_t * [nb_snps];
  size_t k = 0;
  for(size_t i = 0; i < p_A->nrow; i++) {
    if(snps[i]) data[k++] = p_A->data[i];
  }
  paraKin_p X(data, p_A->ncol, p_A->true_ncol, (Ktype) (nb_snps - 1), &p[0], dominance);
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
      Y(i,j) = (double) X.K[k++]; // ou Y(j,i)
    }
  }

  return Y;
}


RcppExport SEXP gg_Kinship_pw(SEXP p_ASEXP, SEXP pSEXP, SEXP snpsSEXP, SEXP domiSEXP, SEXP chunkSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< XPtr<matrix4> >::type p_A(p_ASEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type p(pSEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type snps(snpsSEXP);
    Rcpp::traits::input_parameter< bool >::type domi(domiSEXP);
    Rcpp::traits::input_parameter< int >::type chunk(chunkSEXP);
    __result = Rcpp::wrap(Kinship_pw(p_A, p, snps, domi, chunk));
    return __result;
END_RCPP
}

