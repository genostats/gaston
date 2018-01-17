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

RcppExport SEXP gg_which_duplicated_chr_pos(SEXP ChrSEXP, SEXP PosSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type Chr(ChrSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type Pos(PosSEXP);
    rcpp_result_gen = Rcpp::wrap(which_duplicated_chr_pos(Chr, Pos));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP gg_which_duplicated_chr_pos_alleles(SEXP ChrSEXP, SEXP PosSEXP, SEXP AL1SEXP, SEXP AL2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type Chr(ChrSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type Pos(PosSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type AL1(AL1SEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type AL2(AL2SEXP);
    rcpp_result_gen = Rcpp::wrap(which_duplicated_chr_pos_alleles(Chr, Pos, AL1, AL2));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP gg_which_duplicated_id(SEXP IdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type Id(IdSEXP);
    rcpp_result_gen = Rcpp::wrap(which_duplicated_id(Id));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP gg_which_duplicated_id_chr_pos(SEXP IdSEXP, SEXP ChrSEXP, SEXP PosSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type Id(IdSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type Chr(ChrSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type Pos(PosSEXP);
    rcpp_result_gen = Rcpp::wrap(which_duplicated_id_chr_pos(Id, Chr, Pos));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP gg_which_duplicated_id_chr_pos_alleles(SEXP IdSEXP, SEXP ChrSEXP, SEXP PosSEXP, SEXP AL1SEXP, SEXP AL2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type Id(IdSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type Chr(ChrSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type Pos(PosSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type AL1(AL1SEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type AL2(AL2SEXP);
    rcpp_result_gen = Rcpp::wrap(which_duplicated_id_chr_pos_alleles(Id, Chr, Pos, AL1, AL2));
    return rcpp_result_gen;
END_RCPP
}

