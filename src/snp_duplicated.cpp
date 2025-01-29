#include <Rcpp.h>
#include <iostream>
#include "snp_hash.h"
#include "flip_strand.h"

using namespace Rcpp;

//[[Rcpp::export]]
IntegerVector which_duplicated_chr_pos(IntegerVector Chr, IntegerVector Pos) {
  SNPhash h(Chr, Pos);
  return wrap(h.dup_indices);
}

//[[Rcpp::export]]
IntegerVector which_duplicated_chr_pos_alleles(IntegerVector Chr, IntegerVector Pos, CharacterVector AL1, CharacterVector AL2) {
  SNPhash h(Chr, Pos, AL1, AL2);
  return wrap(h.dup_indices);
}

//[[Rcpp::export]]
IntegerVector which_duplicated_id(CharacterVector Id) {
  SNPhash h(Id);
  return wrap(h.dup_indices);
}

//[[Rcpp::export]]
IntegerVector which_duplicated_id_chr_pos(CharacterVector Id, IntegerVector Chr, IntegerVector Pos) {
  SNPhash h(Id, Chr, Pos);
  return wrap(h.dup_indices);
}

//[[Rcpp::export]]
IntegerVector which_duplicated_id_chr_pos_alleles(CharacterVector Id, IntegerVector Chr, IntegerVector Pos, CharacterVector AL1, CharacterVector AL2) {
  SNPhash h(Id, Chr, Pos, AL1, AL2);
  return wrap(h.dup_indices);
}


