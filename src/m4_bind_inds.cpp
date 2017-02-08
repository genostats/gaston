#include <Rcpp.h>
#include <iostream>
#include "matrix4.h"

using namespace Rcpp;

XPtr<matrix4> bind_inds(List L) {
  int s = L.size();
  if(s < 2) stop("Can't bind less than two matrices!");
  XPtr<matrix4> first = as<XPtr<matrix4> >(L[0]);
  int n = first->ncol;
  int m = first->nrow;
  for(int i = 1; i < s; i++) {
    XPtr<matrix4> nxt = as<XPtr<matrix4> >(L[i]);
    if(m != nxt->nrow) stop("Dimensions mismatch");
    n += nxt->ncol;
  }
  XPtr<matrix4> r(new matrix4(m,n));
  for(int i = 0; i < m; i++) {
    int k = 0;
    for(int j = 0; j < s; j++) {
      XPtr<matrix4> nxt = as<XPtr<matrix4> >(L[j]);
      for(int jj = 0; jj < nxt->ncol; jj++) {
        (*r)(i,k++) = (*nxt)(i,jj); }
    }
  }
  return r;
}

inline char flip_strand(char x) {
  if(x == 'A') return 'T';
  if(x == 'C') return 'G';
  if(x == 'G') return 'C';
  if(x == 'T') return 'A';
  stop("Alleles should be coded A/C/G/T\n");
  return 'X'; // suppress a warning
}

// ci-dessous il faut vérifier en amont que
// * tous les data frames ont des composantes A1 et A2
// * il n'y a que des allèles A C G T
List alleles_recoding(List L) {
  int s = L.size();
  if(s < 2) stop("Can't bind less than two matrices!");
  DataFrame first = as<DataFrame> (L[0]);
  int n = first.nrows();

  int nas = 0, strand = 0, ref = 0;
  
  LogicalMatrix R(s,n);
  for(int i = 0; i < n; i++) { // on parcourt les SNP
    char A1 = as<CharacterVector>(first["A1"])[i][0];  // le premier charactere...
    char A2 = as<CharacterVector>(first["A2"])[i][0];  
    R(0,i) = false;
    for(int j = 1; j < s; j++) {
      DataFrame next = as<DataFrame> (L[j]);
      char B1 = as<CharacterVector>(next["A1"])[i][0]; 
      char B2 = as<CharacterVector>(next["A2"])[i][0];
      // tout ok
      if( A1 == B1 && A2 == B2 ) {
        R(j,i) = false;
        continue; // case closed
      }
      // inversion allèle ref ?
      if( A1 == B2 && A2 == B1 ) {
        R(j,i) = true;
        ref++;
        continue; // case closed
      }
      // remarque : si les deux SNPs sont CG, un des deux cas précédents a réglé la question
      // idem pour deux SNPs AT
      // on passe à une éventuelle inversion de brin
      char C1 = flip_strand(B1);
      char C2 = flip_strand(B2);
      if( A1 == C1 && A2 == C2 ) { // ok
        R(j,i) = false;
        strand++;
        continue;
      }
      if( A1 == C2 && A2 == C1 ) { // + inversion allèle ref
        R(j,i) = true;
        strand++; 
        ref++;
        continue; // case closed
      }
      // si on est arrivés jusque là c'est qu'il y a discordance sans remède
      for(int k = 0; k < s; k++) R(k,i) = NA_LOGICAL;
      nas++;
      break; // on a réglé le cas de toute la colonne
    }
  }
  List LL;
  LL["flip"] = R;
  LL["ref"] = ref;
  LL["strand"] = strand;
  LL["NAs"] = nas;
  return LL;
}

inline uint8_t flip_strand2(uint8_t x) {
  if(x==3) return 3;
  return(2-x);
}

//[[Rcpp::export]]
XPtr<matrix4> bind_inds2(List L, LogicalMatrix flip) {
  int s = L.size();
  if(s < 2) 
    stop("Can't bind less than two matrices!");
  if(flip.nrow() != s)
    stop("Dimensions mismatch");

  XPtr<matrix4> first = as<XPtr<matrix4> >(L[0]);
  int n = first->ncol;
  int m = first->nrow;
  if(flip.ncol() != m)
    stop("Dimensions mismatch");
  for(int i = 1; i < s; i++) {
    XPtr<matrix4> nxt = as<XPtr<matrix4> >(L[i]);
    if(m != nxt->nrow) stop("Dimensions mismatch");
    n += nxt->ncol;
  }
  XPtr<matrix4> r(new matrix4(m,n));
  for(int i = 0; i < m; i++) {
    int k = 0;
    for(int j = 0; j < s; j++) {
      XPtr<matrix4> nxt = as<XPtr<matrix4> >(L[j]);
      for(int jj = 0; jj < nxt->ncol; jj++) {
        if(LogicalVector::is_na(flip(j,i))) {
          (*r)(i,k++) = 3; // NA
        } else if(flip(j,i)) {
          (*r)(i,k++) = flip_strand2( (*nxt)(i,jj) );
        } else {
          (*r)(i,k++) = (*nxt)(i,jj);
        }
      }
    }
  }
  return r;
}









RcppExport SEXP gg_bind_inds(SEXP LSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< List >::type L(LSEXP );
        XPtr<matrix4> __result = bind_inds(L);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}

RcppExport SEXP gg_alleles_recoding(SEXP LSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< List >::type L(LSEXP);
    __result = Rcpp::wrap(alleles_recoding(L));
    return __result;
END_RCPP
}

RcppExport SEXP gg_bind_inds2(SEXP LSEXP, SEXP flipSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< List >::type L(LSEXP);
    Rcpp::traits::input_parameter< LogicalMatrix >::type flip(flipSEXP);
    __result = Rcpp::wrap(bind_inds2(L, flip));
    return __result;
END_RCPP
}

