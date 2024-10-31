#include <Rcpp.h>
#include <iostream>
#include "matrix4.h"
#include "loubar.h"
#include "m4_ld.h"

using namespace Rcpp;

// un matrice de génotypes, les deux vecteurs nécessaires à centrer/réduire
// threshold : le seuil de r² accepté
// pos : les positions des SNPs 
// chr : les chromosomes 
// max_dist : la distance au dessus de laquelle on ne calcule pas le LD
// beg, end : les indices des SNPs à considérer (inclusifs)
// w_ = Logical Vecteur, quels SNPs considérer ? (true = keep)
// [[Rcpp::export]]
LogicalVector ld_thin_random(XPtr<matrix4> pA, NumericVector mu, NumericVector sd, double threshold, 
                      IntegerVector pos, IntegerVector chr, int max_dist, int beg, int end,
                      LogicalVector w_) {
  int n = end - beg + 1;

  LogicalVector w;
  if(w_.length() == 0) {
    w = LogicalVector(n);
    for(int i = 0; i < n; i++) w(i) = true;
  }
  else if(w_.length() == n)
    w = clone(w_);
  else
    stop("Length of which.snps doesn't match\n");

  int i = beg;
  threshold = sqrt(threshold);
  while(i <= end) {
    int j = i + 1;
    int max_pos = pos(i) + max_dist;
    int chr_i = chr(i);
    double mu_i = mu(i);
    double sd_i = sd(i);
    int next_i;
    bool gotnexti = false;
    while(j <= end && pos(j) < max_pos && chr(j) == chr_i) {

      if(!w(j-beg)) { 
        j++; 
        continue; 
      }

      double ld = LD_colxx(*pA, mu_i, mu(j), sd_i*sd(j), i, j);
 
      if(fabs(ld) < threshold) {
        if(!gotnexti) {
          next_i = j;
          gotnexti = true;
        }
        j++;
        continue;
      }
      if(Rf_runif(0,1) < 0.5) {
        // keep left;
        w(j-beg) = false;
        j++;
        continue; // next j...
      }
      // keep right
      w(i-beg) = false; 
      break; // next i
    }
    if(gotnexti) i = next_i; else i = j; 
  }
  return w; 
}
