#include <Rcpp.h>
#include <iostream>
#include <fstream>
#include <ctime>
#include "matrix4.h"
#include "loubar.h"

using namespace Rcpp;


/*** 
  Inversion des allèles
  0 -> 2
  1 -> 1
  2 -> 0
  3 -> NA
***/


uint8_t rec[256] = {
170, 169, 168, 171, 166, 165, 164, 167, 162, 161, 160, 163, 174, 173, 172, 175, 
154, 153, 152, 155, 150, 149, 148, 151, 146, 145, 144, 147, 158, 157, 156, 159, 
138, 137, 136, 139, 134, 133, 132, 135, 130, 129, 128, 131, 142, 141, 140, 143, 
186, 185, 184, 187, 182, 181, 180, 183, 178, 177, 176, 179, 190, 189, 188, 191, 
106, 105, 104, 107, 102, 101, 100, 103,  98,  97,  96,  99, 110, 109, 108, 111, 
 90,  89,  88,  91,  86,  85,  84,  87,  82,  81,  80,  83,  94,  93,  92,  95, 
 74,  73,  72,  75,  70,  69,  68,  71,  66,  65,  64,  67,  78,  77,  76,  79, 
122, 121, 120, 123, 118, 117, 116, 119, 114, 113, 112, 115, 126, 125, 124, 127, 
 42,  41,  40,  43,  38,  37,  36,  39,  34,  33,  32,  35,  46,  45,  44,  47, 
 26,  25,  24,  27,  22,  21,  20,  23,  18,  17,  16,  19,  30,  29,  28,  31, 
 10,   9,   8,  11,   6,   5,   4,   7,   2,   1,   0,   3,  14,  13,  12,  15, 
 58,  57,  56,  59,  54,  53,  52,  55,  50,  49,  48,  51,  62,  61,  60,  63, 
234, 233, 232, 235, 230, 229, 228, 231, 226, 225, 224, 227, 238, 237, 236, 239, 
218, 217, 216, 219, 214, 213, 212, 215, 210, 209, 208, 211, 222, 221, 220, 223, 
202, 201, 200, 203, 198, 197, 196, 199, 194, 193, 192, 195, 206, 205, 204, 207, 
250, 249, 248, 251, 246, 245, 244, 247, 242, 241, 240, 243, 254, 253, 252, 255 };

uint8_t hzna[256] = {
  0,   3,   2,   3,  12,  15,  14,  15,   8,  11,  10,  11,  12,  15,  14,  15, 
 48,  51,  50,  51,  60,  63,  62,  63,  56,  59,  58,  59,  60,  63,  62,  63, 
 32,  35,  34,  35,  44,  47,  46,  47,  40,  43,  42,  43,  44,  47,  46,  47, 
 48,  51,  50,  51,  60,  63,  62,  63,  56,  59,  58,  59,  60,  63,  62,  63, 
192, 195, 194, 195, 204, 207, 206, 207, 200, 203, 202, 203, 204, 207, 206, 207, 
240, 243, 242, 243, 252, 255, 254, 255, 248, 251, 250, 251, 252, 255, 254, 255, 
224, 227, 226, 227, 236, 239, 238, 239, 232, 235, 234, 235, 236, 239, 238, 239, 
240, 243, 242, 243, 252, 255, 254, 255, 248, 251, 250, 251, 252, 255, 254, 255, 
128, 131, 130, 131, 140, 143, 142, 143, 136, 139, 138, 139, 140, 143, 142, 143, 
176, 179, 178, 179, 188, 191, 190, 191, 184, 187, 186, 187, 188, 191, 190, 191, 
160, 163, 162, 163, 172, 175, 174, 175, 168, 171, 170, 171, 172, 175, 174, 175, 
176, 179, 178, 179, 188, 191, 190, 191, 184, 187, 186, 187, 188, 191, 190, 191, 
192, 195, 194, 195, 204, 207, 206, 207, 200, 203, 202, 203, 204, 207, 206, 207, 
240, 243, 242, 243, 252, 255, 254, 255, 248, 251, 250, 251, 252, 255, 254, 255, 
224, 227, 226, 227, 236, 239, 238, 239, 232, 235, 234, 235, 236, 239, 238, 239, 
240, 243, 242, 243, 252, 255, 254, 255, 248, 251, 250, 251, 252, 255, 254, 255 };

void invert_snp_coding(XPtr<matrix4> p_A, size_t snp) {
  if(snp >= p_A->nrow) stop("SNP index out of range");
  uint8_t * d = p_A->data[snp];
  for(size_t j = 0 ; j < p_A->true_ncol; j++) {
    d[j] = rec[ d[j] ];
  }
}


//[[Rcpp::export]]
void snp_hz_to_na(XPtr<matrix4> p_A, size_t snp) {
  if(snp >= p_A->nrow) stop("SNP index out of range");
  uint8_t * d = p_A->data[snp];
  for(size_t j = 0 ; j < p_A->true_ncol; j++) {
    d[j] = hzna[ d[j] ];
  }
}

//[[Rcpp::export]]
void set_snp_to_na(XPtr<matrix4> p_A, size_t snp) {
  if(snp >= p_A->nrow) stop("SNP index out of range");
  uint8_t * d = p_A->data[snp];
  for(size_t j = 0 ; j < p_A->true_ncol; j++) {
    d[j] = 255;
  }
}



RcppExport SEXP gg_invert_snp_coding(SEXP p_ASEXP, SEXP snpSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< XPtr<matrix4> >::type p_A(p_ASEXP);
    Rcpp::traits::input_parameter< size_t >::type snp(snpSEXP);
    invert_snp_coding(p_A, snp);
    return R_NilValue;
END_RCPP
}

RcppExport SEXP gg_snp_hz_to_na(SEXP p_ASEXP, SEXP snpSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< XPtr<matrix4> >::type p_A(p_ASEXP);
    Rcpp::traits::input_parameter< size_t >::type snp(snpSEXP);
    snp_hz_to_na(p_A, snp);
    return R_NilValue;
END_RCPP
}

RcppExport SEXP gg_set_snp_to_na(SEXP p_ASEXP, SEXP snpSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< XPtr<matrix4> >::type p_A(p_ASEXP);
    Rcpp::traits::input_parameter< size_t >::type snp(snpSEXP);
    set_snp_to_na(p_A, snp);
    return R_NilValue;
END_RCPP
}

