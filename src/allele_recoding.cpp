#include <Rcpp.h>
#include <iostream>
#include "matrix4.h"
#include "flip_strand.h"

using namespace Rcpp;

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
    CharacterVector A1 = as<CharacterVector>(first["A1"]);
    CharacterVector A2 = as<CharacterVector>(first["A2"]);
    const char * a1 = CHAR(STRING_ELT(A1,i));
    const char * a2 = CHAR(STRING_ELT(A2,i));
    R(0,i) = false;
    for(int j = 1; j < s; j++) {
      DataFrame next = as<DataFrame> (L[j]);
      CharacterVector B1 = as<CharacterVector>(next["A1"]);
      CharacterVector B2 = as<CharacterVector>(next["A2"]);
      const char * b1 = CHAR(STRING_ELT(B1,i));
      const char * b2 = CHAR(STRING_ELT(B2,i));
      // tout ok
      if( !std::strcmp(a1, b1) && !std::strcmp(a2, b2) ) {
        R(j,i) = false;
        continue; // case closed
      }
      // inversion allèle ref ?
      if( !std::strcmp(a1, b2) && !std::strcmp(a2, b1) ) {
        R(j,i) = true;
        ref++;
        continue; // case closed
      }
      // remarque : si les deux SNPs sont CG, un des deux cas précédents a réglé la question
      // idem pour deux SNPs AT
      // on passe à une éventuelle inversion de brin
      std::string c1 = flip_strand(b1);
      std::string c2 = flip_strand(b2);
      if( c1 == a1 && c2 == a2 ) { // ok
        R(j,i) = false;
        strand++;
        continue;
      }
      if( c1 == a2 && c2 == a1 ) { // + inversion allèle ref
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
  LL["swap_reference"] = ref;
  LL["flip_strand"] = strand;
  LL["NAs"] = nas;
  return LL;
}

// -----------------------------------------------------------------

RcppExport SEXP gg_alleles_recoding(SEXP LSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< List >::type L(LSEXP);
    __result = Rcpp::wrap(alleles_recoding(L));
    return __result;
END_RCPP
}
