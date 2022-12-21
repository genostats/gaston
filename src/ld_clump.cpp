#include <Rcpp.h>
#include "matrix4.h"
#include "m4_ld.h"

using namespace Rcpp;

// un matrice de génotypes, les deux vecteurs nécessaires à centrer/réduire
// threshold : le seuil de r² accepté
// pos : les positions des SNPs 
// chr : les chromosomes 
// max_dist : la distance au dessus de laquelle on ne calcule pas le LD
// order : l'ordre dans lequel considérer les SNPs = order(p) - 1
// 
// !!! SUPPOSE QUE LES SNPs SONT BIEN DANS L'ORDRE DU GENOME !!!
//
// [[Rcpp::export]]
IntegerVector ld_clump(XPtr<matrix4> pA, NumericVector mu, NumericVector sd, double threshold,
                       IntegerVector pos, IntegerVector chr, int max_dist, IntegerVector order) {

  IntegerVector Index(pA->nrow, -1);
  threshold = sqrt(threshold); //!!!!!!!!!!!!! on passe de r² à |r| !!!!!!!!!!!!!!!!

  for(int i : order) {
    if( Index[i] >= 0 ) // déjà indexé
      continue;
    // sinon ce SNP devient un index, et on cherche tous ceux qui sont 
    // attrapés par lui
    Index[i] = i;
    int c = chr[i];
    int min_pos = pos[i] - max_dist;
    int max_pos = pos[i] + max_dist;
    double mu_i = mu[i];
    double sd_i = sd[i];
    // on balaie le chr vers la gauche
    int j = i-1;
    while(chr[j] == c && pos[j] > min_pos) {
      if(Index[j] < 0) { // SNP pas encore indexé
        double ld = LD_colxx(*pA, mu_i, mu[j], sd_i*sd[j], i, j);
        if(fabs(ld) >= threshold) // si il est en LD avec SNP i
          Index[j] = i;          // indexé par i
      }
      j--;
    }
    // puis vers la droite
    j = i+1;
    while(chr[j] == c && pos[j] < max_pos) {
      if(Index[j] < 0) {
        double ld = LD_colxx(*pA, mu_i, mu[j], sd_i*sd[j], i, j);
        if(fabs(ld) >= threshold)
          Index[j] = i; 
      }
      j++;
    }
  }
  return Index;
}

RcppExport SEXP ld_clump(SEXP pASEXP, SEXP muSEXP, SEXP sdSEXP, SEXP thresholdSEXP, SEXP posSEXP, SEXP chrSEXP, SEXP max_distSEXP, SEXP orderSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< XPtr<matrix4> >::type pA(pASEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sd(sdSEXP);
    Rcpp::traits::input_parameter< double >::type threshold(thresholdSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type pos(posSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type chr(chrSEXP);
    Rcpp::traits::input_parameter< int >::type max_dist(max_distSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type order(orderSEXP);
    rcpp_result_gen = Rcpp::wrap(ld_clump(pA, mu, sd, threshold, pos, chr, max_dist, order));
    return rcpp_result_gen;
END_RCPP
}


