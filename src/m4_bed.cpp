#include <Rcpp.h>
#include <iostream>
#include <fstream>
#include "matrix4.h"

using namespace Rcpp;


/*** 
  Dans le format bed : 00 et 11 pour les deux homoz 
                       10 pour hétéroz
                       01 pour manquant

  Donc il faut une conversion 0 -> 0 ou 00 -> 00
                              1 -> 3    01 -> 11
                              2 -> 1    10 -> 01
                              3 -> 2    11 -> 10 

  on peut le faire simultanément sur les 8 génotypes contenus
  dans un uint8_t
***/

uint8_t bedc[256] = {
  0,   3,   1,   2,  12,  15,  13,  14,   4,   7,   5,   6,   8,  11,   9,  10, 
 48,  51,  49,  50,  60,  63,  61,  62,  52,  55,  53,  54,  56,  59,  57,  58, 
 16,  19,  17,  18,  28,  31,  29,  30,  20,  23,  21,  22,  24,  27,  25,  26, 
 32,  35,  33,  34,  44,  47,  45,  46,  36,  39,  37,  38,  40,  43,  41,  42, 
192, 195, 193, 194, 204, 207, 205, 206, 196, 199, 197, 198, 200, 203, 201, 202, 
240, 243, 241, 242, 252, 255, 253, 254, 244, 247, 245, 246, 248, 251, 249, 250, 
208, 211, 209, 210, 220, 223, 221, 222, 212, 215, 213, 214, 216, 219, 217, 218, 
224, 227, 225, 226, 236, 239, 237, 238, 228, 231, 229, 230, 232, 235, 233, 234, 
 64,  67,  65,  66,  76,  79,  77,  78,  68,  71,  69,  70,  72,  75,  73,  74,  
112, 115, 113, 114, 124, 127, 125, 126, 116, 119, 117, 118, 120, 123, 121, 122, 
 80,  83,  81,  82,  92,  95,  93,  94,  84,  87,  85,  86,  88,  91,  89,  90,  
 96,  99,  97,  98, 108, 111, 109, 110, 100, 103, 101, 102, 104, 107, 105, 106, 
128, 131, 129, 130, 140, 143, 141, 142, 132, 135, 133, 134, 136, 139, 137, 138, 
176, 179, 177, 178, 188, 191, 189, 190, 180, 183, 181, 182, 184, 187, 185, 186, 
144, 147, 145, 146, 156, 159, 157, 158, 148, 151, 149, 150, 152, 155, 153, 154, 
160, 163, 161, 162, 172, 175, 173, 174, 164, 167, 165, 166, 168, 171, 169, 170};


uint8_t tobed[256] = {
  0,   2,   3,   1,   8,  10,  11,   9,  12,  14,  15,  13,   4,   6,   7,   5,
 32,  34,  35,  33,  40,  42,  43,  41,  44,  46,  47,  45,  36,  38,  39,  37,
 48,  50,  51,  49,  56,  58,  59,  57,  60,  62,  63,  61,  52,  54,  55,  53,
 16,  18,  19,  17,  24,  26,  27,  25,  28,  30,  31,  29,  20,  22,  23,  21,
128, 130, 131, 129, 136, 138, 139, 137, 140, 142, 143, 141, 132, 134, 135, 133,
160, 162, 163, 161, 168, 170, 171, 169, 172, 174, 175, 173, 164, 166, 167, 165,
176, 178, 179, 177, 184, 186, 187, 185, 188, 190, 191, 189, 180, 182, 183, 181,
144, 146, 147, 145, 152, 154, 155, 153, 156, 158, 159, 157, 148, 150, 151, 149,
192, 194, 195, 193, 200, 202, 203, 201, 204, 206, 207, 205, 196, 198, 199, 197,
224, 226, 227, 225, 232, 234, 235, 233, 236, 238, 239, 237, 228, 230, 231, 229,
240, 242, 243, 241, 248, 250, 251, 249, 252, 254, 255, 253, 244, 246, 247, 245,
208, 210, 211, 209, 216, 218, 219, 217, 220, 222, 223, 221, 212, 214, 215, 213,
 64,  66,  67,  65,  72,  74,  75,  73,  76,  78,  79,  77,  68,  70,  71,  69,
 96,  98,  99,  97, 104, 106, 107, 105, 108, 110, 111, 109, 100, 102, 103, 101,
112, 114, 115, 113, 120, 122, 123, 121, 124, 126, 127, 125, 116, 118, 119, 117,
 80,  82,  83,  81,  88,  90,  91,  89,  92,  94,  95,  93,  84,  86,  87,  85};

// bed magic numbers : 108 27 1

XPtr<matrix4> read_bed_file(CharacterVector filename, int n_ind, int n_snp) {
  std::ifstream file(filename[0], std::ifstream::binary);
  if(!file.is_open()) {
    Rf_error("Cannot open file");
  }
  uint8_t m1, m2, m3;
  file.read(reinterpret_cast<char *>(&m1), 1);
  file.read(reinterpret_cast<char *>(&m2), 1);
  file.read(reinterpret_cast<char *>(&m3), 1);
  if(m1 != 108 || m2 != 27) {
    Rf_error("Not a bed file");
  }
  if(m3 != 1) {
    Rf_error("Not a bed file in SNP major mode");
  }
  XPtr<matrix4> p_A(new matrix4(n_snp, n_ind));
  uint8_t b;
  uint8_t bordermask;
  switch(4*p_A->true_ncol - n_ind) {
    case 0:
      bordermask = 0;
      break;
    case 1:
      bordermask = 192;
      break;
    case 2: 
      bordermask = 240;
      break;
    case 3:
      bordermask = 252;
      break;
    default:
      Rf_error("Some shit hit the fan very hard");
  }
  for(int i = 0; i < n_snp; i++) {
    for(int j = 0; j < p_A->true_ncol; j++) {
      file.read(reinterpret_cast<char *> (&b),1);
      p_A->data[i][j] = bedc[ (int) b ];
    }
    // être sûr d'être bordé de 11 -> NA
    p_A->data[i][p_A->true_ncol - 1] |= bordermask;
  }
  
  file.close();
  return p_A;
}

// [[Rcpp::export]]
void write_bed_file(XPtr<matrix4> p_A, CharacterVector filename) {
  std::ofstream file(filename[0], std::ofstream::binary);
  if(!file.is_open()) {
    Rf_error("Cannot open file");
  }

  uint8_t magic[3] = {108, 27, 1};
  file.write(reinterpret_cast<char *> (magic), 3);

  int n_snp = p_A->nrow;
  uint8_t b;
  for(int i = 0; i < n_snp; i++) {
    for(int j = 0; j < p_A->true_ncol; j++) {
      b = tobed[ (int) p_A->data[i][j] ];
      file.write( reinterpret_cast<char *> (&b),1);
    }
  }

  file.close();
}





RcppExport SEXP gg_read_bed_file(SEXP filenameSEXP, SEXP n_indSEXP, SEXP n_snpSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< CharacterVector >::type filename(filenameSEXP );
        Rcpp::traits::input_parameter< int >::type n_ind(n_indSEXP );
        Rcpp::traits::input_parameter< int >::type n_snp(n_snpSEXP );
        XPtr<matrix4> __result = read_bed_file(filename, n_ind, n_snp);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}



RcppExport SEXP gg_write_bed_file(SEXP p_ASEXP, SEXP filenameSEXP) {
BEGIN_RCPP
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< XPtr<matrix4> >::type p_A(p_ASEXP );
        Rcpp::traits::input_parameter< CharacterVector >::type filename(filenameSEXP );
        write_bed_file(p_A, filename);
    }
    return R_NilValue;
END_RCPP
}

