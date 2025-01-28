#include <Rcpp.h>
#include <iostream>
#include "matrix4.h"
#include "ROH.h"

using namespace Rcpp;

// pos is either positions in kb or distance from beginning in cM...
// [[Rcpp::export]]
List ROHs(XPtr<matrix4> pA, IntegerVector chr, NumericVector pos, int beg, int end, double minROHLength, double minDistHet, double maxGapLength, bool NAsAreHet) {
  size_t ncol = pA->ncol;
  std::vector<ROH<ROHsegments>> R(ncol);

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

  // wrapping the result
  List L(ncol);
  for(unsigned int i = 0; i < ncol; i++) {
    DataFrame D = DataFrame::create(Named("i.beg") = wrap(R[i].summary.begIndex), Named("i.end") = wrap(R[i].summary.endIndex), 
        Named("pos.beg") = wrap(R[i].summary.begPos), Named("pos.end") = wrap(R[i].summary.endPos), 
        Named("nb.het") = wrap(R[i].summary.nbHet));
    L[i] = D;
  }
  return L;
}

