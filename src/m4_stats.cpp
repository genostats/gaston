#include <Rcpp.h>
#include <iostream>
#include <ctime>
#include "matrix4.h"

using namespace Rcpp;

uint8_t N0[256] = {
4, 3, 3, 3, 3, 2, 2, 2, 3, 2, 2, 2, 3, 2, 2, 2, 
3, 2, 2, 2, 2, 1, 1, 1, 2, 1, 1, 1, 2, 1, 1, 1, 
3, 2, 2, 2, 2, 1, 1, 1, 2, 1, 1, 1, 2, 1, 1, 1, 
3, 2, 2, 2, 2, 1, 1, 1, 2, 1, 1, 1, 2, 1, 1, 1, 
3, 2, 2, 2, 2, 1, 1, 1, 2, 1, 1, 1, 2, 1, 1, 1, 
2, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 
2, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 
2, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 
3, 2, 2, 2, 2, 1, 1, 1, 2, 1, 1, 1, 2, 1, 1, 1, 
2, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 
2, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 
2, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 
3, 2, 2, 2, 2, 1, 1, 1, 2, 1, 1, 1, 2, 1, 1, 1, 
2, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 
2, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 
2, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0};

uint8_t N1[256] = {
0, 1, 0, 0, 1, 2, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 
1, 2, 1, 1, 2, 3, 2, 2, 1, 2, 1, 1, 1, 2, 1, 1, 
0, 1, 0, 0, 1, 2, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 
0, 1, 0, 0, 1, 2, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 
1, 2, 1, 1, 2, 3, 2, 2, 1, 2, 1, 1, 1, 2, 1, 1, 
2, 3, 2, 2, 3, 4, 3, 3, 2, 3, 2, 2, 2, 3, 2, 2, 
1, 2, 1, 1, 2, 3, 2, 2, 1, 2, 1, 1, 1, 2, 1, 1, 
1, 2, 1, 1, 2, 3, 2, 2, 1, 2, 1, 1, 1, 2, 1, 1, 
0, 1, 0, 0, 1, 2, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 
1, 2, 1, 1, 2, 3, 2, 2, 1, 2, 1, 1, 1, 2, 1, 1, 
0, 1, 0, 0, 1, 2, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 
0, 1, 0, 0, 1, 2, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 
0, 1, 0, 0, 1, 2, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 
1, 2, 1, 1, 2, 3, 2, 2, 1, 2, 1, 1, 1, 2, 1, 1, 
0, 1, 0, 0, 1, 2, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 
0, 1, 0, 0, 1, 2, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0};


// TOUT SUPpOSE QUE LA MATRICE EST BORDEE DE 3/na

// Stats SNP : NO N1 N2 tous individus confondus / seulement chez les femmes pour le chr X
// Stats Individus : callrate et hz autosomal, chr X, chr Y, mt séparés
List geno_stats(matrix4 & A, LogicalVector chr_x, LogicalVector chr_y, LogicalVector chr_mt, LogicalVector sexf) {
  // création vecteur de masques sex
  uint8_t * sex_ = new uint8_t[A.true_ncol];
  std::fill(sex_, sex_+A.true_ncol, 0);  
  int nbm = 0;
  for(size_t j = 0; j < A.true_ncol; j++) {
    if(!sexf(4*j)) { sex_[j] |= 3; nbm++; }
    if(4*j+1 < A.ncol && !sexf(4*j+1)) { sex_[j] |= 12;  nbm++; }
    if(4*j+2 < A.ncol && !sexf(4*j+2)) { sex_[j] |= 48;  nbm++; }
    if(4*j+3 < A.ncol && !sexf(4*j+3)) { sex_[j] |= 192; nbm++; }
  }

  // Stats SNP
  IntegerMatrix SN(4,A.nrow);
  IntegerMatrix SNf(4,A.nrow); // à remplir seulement pour chr X

  // Stats individus
  IntegerMatrix sIN(4,   A.true_ncol*4);  // autosomes
  IntegerMatrix sINx(4,  A.true_ncol*4);  // chr x
  IntegerMatrix sINy(4,  A.true_ncol*4);  // chr y 
  IntegerMatrix sINmt(4, A.true_ncol*4);  // chr mt

  int * pIN   = sIN.begin();
  int * pINx  = sINx.begin();
  int * pINy  = sINy.begin();
  int * pINmt = sINmt.begin();

  for(size_t i = 0; i < A.nrow; i++) {
    // Quelle stat individus à mettre jour, selon chromosome
    int * pIN_; 
    bool xy = false;
    if(chr_x(i)) { 
      pIN_ = pINx;
      xy = true;
    } else if(chr_y(i)) {
      pIN_ = pINy;
      xy = true;
    } else if(chr_mt(i))
      pIN_ = pINmt;
    else
      pIN_ = pIN;

    for(size_t j = 0; j < A.true_ncol; j++) {
      uint8_t d = A.data[i][j];
      SN(0,i) += N0[d];
      SN(3,i) += N0[255-d];
      SN(1,i) += N1[d];
      SN(2,i) += N1[255-d];

      // stat Individus
      pIN_[16*j + ((int) d&3)]++;
      pIN_[16*j + 4 + (int) ((d>>2)&3)]++;
      pIN_[16*j + 8 + (int) ((d>>4)&3)]++;
      pIN_[16*j + 12 + (int) ((d>>6)&3)]++;

      if(xy) { // Stats SNPs en ne prenant en compte que les femmes
        d |= sex_[j];
        SNf(0,i) += N0[d];
        SNf(3,i) += N0[255-d];
        SNf(1,i) += N1[d];
        SNf(2,i) += N1[255-d];
      }
    }
    // Des NAs en trop (la bordure)
    SN(3,i) -= (4*A.true_ncol - A.ncol);
    // Pour les femmes, 
    // le masque met les hommes à NA et ne touche pas à la bordure
    SNf(3,i) -= (4*A.true_ncol - A.ncol) + nbm;
  }
  
  sIN   = sIN(_, Range(0,A.ncol-1));
  sINx  = sINx(_, Range(0,A.ncol-1));
  sINy  = sINy(_, Range(0,A.ncol-1));
  sINmt = sINmt(_, Range(0,A.ncol-1));

  List L;

  L["snps"] = DataFrame::create(Named("N0") = SN(0,_), Named("N1")  = SN(1,_), 
                                Named("N2") = SN(2,_), Named("NAs") = SN(3,_),
                                Named("N0.f") = SNf(0,_), Named("N1.f") = SNf(1,_),
                                Named("N2.f") = SNf(2,_), Named("NAs.f") = SNf(3,_));

  L["inds"] = DataFrame::create(Named("N0") = sIN(0,_), Named("N1")  = sIN(1,_), 
                                Named("N2") = sIN(2,_), Named("NAs") = sIN(3,_),
                                Named("N0.x") = sINx(0,_), Named("N1.x")  = sINx(1,_), 
                                Named("N2.x") = sINx(2,_), Named("NAs.x") = sINx(3,_),
                                Named("N0.y") = sINy(0,_), Named("N1.y")  = sINy(1,_), 
                                Named("N2.y") = sINy(2,_), Named("NAs.y") = sINy(3,_),
                                Named("N0.mt") = sINmt(0,_), Named("N1.mt")  = sINmt(1,_), 
                                Named("N2.mt") = sINmt(2,_), Named("NAs.mt") = sINmt(3,_) );
  delete [] sex_;
  return L;
}

List geno_stats_snps(matrix4 & A, LogicalVector chr_xy, LogicalVector sexf) {

  // création vecteur de masques sex
  uint8_t * sex_ = new uint8_t[A.true_ncol];
  std::fill(sex_, sex_+A.true_ncol, 0);  
  int nbm = 0;
  for(size_t j = 0; j < A.true_ncol; j++) {
    if(!sexf(4*j)) { sex_[j] |= 3; nbm++; }
    if(4*j+1 < A.ncol && !sexf(4*j+1)) { sex_[j] |= 12;  nbm++; }
    if(4*j+2 < A.ncol && !sexf(4*j+2)) { sex_[j] |= 48;  nbm++; }
    if(4*j+3 < A.ncol && !sexf(4*j+3)) { sex_[j] |= 192; nbm++; }
  }

  // Stats SNP
  IntegerMatrix SN(4,A.nrow);
  IntegerMatrix SNf(4,A.nrow); // à remplir seulement pour chr X

  for(size_t i = 0; i < A.nrow; i++) {
    bool xy = chr_xy(i);
    for(size_t j = 0; j < A.true_ncol; j++) {
      uint8_t d = A.data[i][j];
      SN(0,i) += N0[d];
      SN(3,i) += N0[255-d];
      SN(1,i) += N1[d];
      SN(2,i) += N1[255-d];

      if(xy) { // Stats SNPs en ne prenant en compte que les femmes
        d |= sex_[j];
        SNf(0,i) += N0[d];
        SNf(3,i) += N0[255-d];
        SNf(1,i) += N1[d];
        SNf(2,i) += N1[255-d];
      }
    }
    // Des NAs en trop (la bordure)
    SN(3,i) -= (4*A.true_ncol - A.ncol);
    // Pour les femmes la bordure + le nbre d'hommes
    SNf(3,i) -= (4*A.true_ncol - A.ncol) + nbm;
  }
  List L;

  L["snps"] = DataFrame::create(Named("N0") = SN(0,_), Named("N1")  = SN(1,_), 
                                Named("N2") = SN(2,_), Named("NAs") = SN(3,_),
                                Named("N0.f") = SNf(0,_), Named("N1.f") = SNf(1,_),
                                Named("N2.f") = SNf(2,_), Named("NAs.f") = SNf(3,_));

  delete [] sex_;
  return L;
}

List geno_stats_inds(matrix4 & A, LogicalVector chr_x, LogicalVector chr_y, LogicalVector chr_mt) {
  // Stats individus
  IntegerMatrix sIN(4,   A.true_ncol*4);  // autosomes
  IntegerMatrix sINx(4,  A.true_ncol*4);  // chr x
  IntegerMatrix sINy(4,  A.true_ncol*4);  // chr y 
  IntegerMatrix sINmt(4, A.true_ncol*4);  // chr mt

  int * pIN   = sIN.begin();
  int * pINx  = sINx.begin();
  int * pINy  = sINy.begin();
  int * pINmt = sINmt.begin();

  for(size_t i = 0; i < A.nrow; i++) {
    // Quelle stat individus à mettre jour, selon chromosome
    int * pIN_; 
    if(chr_x(i)) { 
      pIN_ = pINx;
    } else if(chr_y(i))
      pIN_ = pINy;
    else if(chr_mt(i))
      pIN_ = pINmt;
    else
      pIN_ = pIN;

    for(size_t j = 0; j < A.true_ncol; j++) {
      uint8_t d = A.data[i][j];
      // stat Individus
      pIN_[16*j + ((int) d&3)]++;
      pIN_[16*j + 4 + (int) ((d>>2)&3)]++;
      pIN_[16*j + 8 + (int) ((d>>4)&3)]++;
      pIN_[16*j + 12 + (int) ((d>>6)&3)]++;
    }
  }
  
  sIN   = sIN(_, Range(0,A.ncol-1));
  sINx  = sINx(_, Range(0,A.ncol-1));
  sINy  = sINy(_, Range(0,A.ncol-1));
  sINmt = sINmt(_, Range(0,A.ncol-1));

  List L;

  L["inds"] = DataFrame::create(Named("N0") = sIN(0,_), Named("N1")  = sIN(1,_), 
                                Named("N2") = sIN(2,_), Named("NAs") = sIN(3,_),
                                Named("N0.x") = sINx(0,_), Named("N1.x")  = sINx(1,_), 
                                Named("N2.x") = sINx(2,_), Named("NAs.x") = sINx(3,_),
                                Named("N0.y") = sINy(0,_), Named("N1.y")  = sINy(1,_), 
                                Named("N2.y") = sINy(2,_), Named("NAs.y") = sINy(3,_),
                                Named("N0.mt") = sINmt(0,_), Named("N1.mt")  = sINmt(1,_), 
                                Named("N2.mt") = sINmt(2,_), Named("NAs.mt") = sINmt(3,_) );
  return L;
}


//[[Rcpp::export]]
List geno_stats(XPtr<matrix4> p_A,  LogicalVector chr_x, LogicalVector chr_y, LogicalVector chr_mt, LogicalVector sexf) {
  return geno_stats(*p_A, chr_x, chr_y, chr_mt, sexf);
}

//[[Rcpp::export]]
List geno_stats_snps(XPtr<matrix4> p_A,  LogicalVector chr_xy, LogicalVector sexf) {
  return geno_stats_snps(*p_A, chr_xy, sexf);
}

//[[Rcpp::export]]
List geno_stats_inds(XPtr<matrix4> p_A,  LogicalVector chr_x, LogicalVector chr_y, LogicalVector chr_mt) {
  return geno_stats_inds(*p_A, chr_x, chr_y, chr_mt);
}



RcppExport SEXP gg_geno_stats(SEXP p_ASEXP, SEXP chr_xSEXP, SEXP chr_ySEXP, SEXP chr_mtSEXP, SEXP sexfSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< XPtr<matrix4> >::type p_A(p_ASEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type chr_x(chr_xSEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type chr_y(chr_ySEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type chr_mt(chr_mtSEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type sexf(sexfSEXP);
    __result = Rcpp::wrap(geno_stats(p_A, chr_x, chr_y, chr_mt, sexf));
    return __result;
END_RCPP
}
// geno_stats_snp
RcppExport SEXP gg_geno_stats_snps(SEXP p_ASEXP, SEXP chr_xySEXP, SEXP sexfSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< XPtr<matrix4> >::type p_A(p_ASEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type chr_xy(chr_xySEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type sexf(sexfSEXP);
    __result = Rcpp::wrap(geno_stats_snps(p_A, chr_xy, sexf));
    return __result;
END_RCPP
}
// geno_stats_ind
RcppExport SEXP gg_geno_stats_inds(SEXP p_ASEXP, SEXP chr_xSEXP, SEXP chr_ySEXP, SEXP chr_mtSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< XPtr<matrix4> >::type p_A(p_ASEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type chr_x(chr_xSEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type chr_y(chr_ySEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type chr_mt(chr_mtSEXP);
    __result = Rcpp::wrap(geno_stats_inds(p_A, chr_x, chr_y, chr_mt));
    return __result;
END_RCPP
}

