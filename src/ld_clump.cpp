#include <Rcpp.h>
#include "matrix4.h"
#include "m4_ld.h"

#include <iostream>

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
  int M = chr.size();

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
    while(j > 0 && chr[j] == c && pos[j] > min_pos) {
      if(Index[j] < 0) { // SNP pas encore indexé
        double ld = LD_colxx(*pA, mu_i, mu[j], sd_i*sd[j], i, j);
        if(fabs(ld) >= threshold) // si il est en LD avec SNP i
          Index[j] = i;          // indexé par i
      }
      j--;
    }
    // puis vers la droite
    j = i+1;
    while(j < M && chr[j] == c && pos[j] < max_pos) {
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

