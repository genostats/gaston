#include <Rcpp.h>
#include <iostream>
#include "matrix4.h"
#include "ROH.h"

using namespace Rcpp;

// pos is either positions in kb or distance from beginning in cM...
// The function will send back a list with two, the total length of ROHs in number SNPs and on the 'pos' scale
// This function is very similar to ROH except for the ROH<ROHlength> at the beginning and the final wrapping of the result
// [[Rcpp::export]]
List ROHlen(XPtr<matrix4> pA, IntegerVector chr, NumericVector pos, int beg, int end, double minROHLength, double minDistHet, double maxGapLength, bool NAsAreHet) {
  size_t ncol = pA->ncol;
  std::vector<ROH<ROHlength>> R(ncol);

  unsigned int current_chr = chr[beg];
  double current_pos = pos[beg];
  for(unsigned int i = beg; i < end; i++) {
    if(chr[i] != current_chr || pos[i] > current_pos + maxGapLength) {
      for(unsigned k = 0; k < ncol; k++) R[k].endChromosome(minROHLength, minDistHet);
      current_chr = chr[i];
    }
    current_pos = pos[i];
    unsigned int k = 0;
    for(unsigned int j = 0; j < pA->true_ncol; j++) {
      uint8_t x = pA->data[i][j];
      for(unsigned int ss = 0; ss < 4 && 4*j + ss < ncol; ss++) {
        uint8_t g = (x&3);
        x >>= 2;
        // le (i+1) ci cessous pour obtenir des résultats avec "R index"
        R[k].update(i+1, pos[i], g, minROHLength, minDistHet, NAsAreHet, false);
        k++;
      }
    }
  }
  // il faut cloturer à la fin du dernier chromosome considéré
  for(unsigned k = 0; k < ncol; k++) R[k].endChromosome(minROHLength, minDistHet);

  // Wrapping the result
  NumericVector nbSNPs(ncol);
  NumericVector length(ncol);
  for(unsigned int i = 0; i < ncol; i++) {
    nbSNPs[i] = R[i].summary.nbSNPs;
    length[i] = R[i].summary.length;
  }
  List L;
  L["nbSNPs"] = nbSNPs;
  L["length"] = length;
  return L;
}

