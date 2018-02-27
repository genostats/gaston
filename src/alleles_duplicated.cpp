#include <Rcpp.h>
#include <iostream>
#include "flip_strand.h"

using namespace Rcpp;

// ci-dessous il faut vérifier en amont que
// * tous les data frames ont des composantes A1 et A2
// * il n'y a que des allèles A C G T
// [HP] pas la peine, si il y a autre chose que des lettres ACGT
// [HP] flip_strand les laisse invariantes, eg flip_strand("ACZ") = "TGZ"
// [HP] donc ça fait pas grand mal
// [HP] D = indices des dupliqués
List alleles_duplicated(DataFrame snps, NumericVector D) {
  int n = snps.nrows();

  int nas = 0, strand = 0, ref = 0;
  
  LogicalVector K(n, true);
  LogicalVector R(n);
  CharacterVector A1 = as<CharacterVector>(snps["A1"]);
  CharacterVector A2 = as<CharacterVector>(snps["A2"]);
  
  for(int i = 0; i < n; i++) { // on parcourt les SNP
  
    if ( R_IsNA(D(i)) ) {
      K(i) = true;
      R(i) = false;
      continue; // pas de dupliqué
    }
    if ( K(i) == false ) continue; //déjà traité
  
    K(i) = true;
    R(i) = false;
  
    const char * a1 = CHAR(STRING_ELT(A1,i));
    const char * a2 = CHAR(STRING_ELT(A2,i));
  
    for(int j = i+1; j < n; j++) { // on cherche les duplica
      if( D(i) != D(j) ) {
        continue;
      }
      const char * b1 = CHAR(STRING_ELT(A1,j));
      const char * b2 = CHAR(STRING_ELT(A2,j));
      // tout ok ?
      if( !std::strcmp(a1, b1) && !std::strcmp(a2, b2) ) {
        K(j) = false;
        R(j) = false;
        continue; // case closed
      }
      // inversion allèle ref ?
      if( !std::strcmp(a1, b2) && !std::strcmp(a2, b1) ) {
        K(j) = false;
        R(j) = true;
        ref++;
        continue; // case closed
      }
      // inversion de brin ?
      std::string c1 = flip_strand(b1);
      std::string c2 = flip_strand(b2);
      if( c1 == a1 && c2 == a2 ) { // ok
        K(j) = false;
        R(j) = false;
        strand++;
        continue;
      }
      if( c1 == a2 && c2 == a1 ) { // + inversion allèle ref
        K(j) = false;
        R(j) = true;
        strand++; 
        ref++;
        continue;
      }
      // un des cas monomorphes
      // [HP] Cas pour SNP avec allèles A/0 vs A/C par exemple... 
      // [HP] (cas rencontré quand plink a lu un fichier ped où il y a 
      // [HP]  A A sur toute la colonnes -> création d'allèles A 0)
      if( !std::strcmp(a1, "0") ) {
        if( !std::strcmp(a2, b2) ) {
          K(i) = false;
          K(j) = true;
          R(j) = false;
          // Il vaut mieux garder le second qui repassera dans la boucle
          // pour un potentiel autre duplica
           break;
        }      
        if( !std::strcmp(a2, b1) ) {
          K(i) = false;
          R(i) = true;
          K(j) = true;
          R(j) = false;
          break;
        }
        std::string c2 = flip_strand(b2);
        if( a2 == c2 ) {
          R(i) = false;
          K(j) = true;
          R(j) = false;
          break;
        }
        std::string c1 = flip_strand(b1);
        if( a2==c1 ) {
          K(i) = false;
          R(i) = true;
          K(j) = true;
          R(j) = false;
          break;
        }
      }
      if( !std::strcmp(a2, "0") ) {
        if( !std::strcmp(a1, b1) ) {
          R(i) = false;
          K(j) = true;
          R(j) = false;
          break;
        }      
        if( !std::strcmp(a1, b2) ) {
          K(i) = false;
          R(i) = true;
          K(j) = true;
          R(j) = false;
          break;
        }
        std::string c1 = flip_strand(b1);
        if( a1==c1 ) {
          R(i) = false;
          K(j) = true;
          R(j) = false;
          break;
        }
        std::string c2 = flip_strand(b2);
        if( a1 == c2 ) {
          K(i) = false;
          R(i) = true;
          K(j) = true;
          R(j) = false;
          break;
        }
      }
      if( !std::strcmp(b1, "0") ) {
        if( !std::strcmp(a2, b2) ) {
          K(j) = false;
          R(j) = false;
          continue;
        }      
        if( !std::strcmp(a1, b2) ) {
          K(j) = false;
          R(j) = true;
          continue;
        }
        std::string c2 = flip_strand(a2);
        if( c2==b2 ) {
          K(j) = false;
          R(j) = false;
         continue;
        }
        std::string c1 = flip_strand(a1);
        if( c1==b2 ) {
          K(j) = false;
          R(j) = true;
         continue;
        }
      }
      if( !std::strcmp(b2, "0") ) {
        if( !std::strcmp(a1, b1) ) {
          K(j) = false;
          R(j) = false;
          continue;
        }      
        if( !std::strcmp(a2, b1) ) {
          K(j) = false;
          R(j) = true;
          continue;
        }
        std::string c1 = flip_strand(a1);
        if( c1==b1 ) {
          K(j) = false;
          R(j) = false;
          continue;
        }
        std::string c2 = flip_strand(a2);
        if( c2==b1 ) {
          K(j) = false;
          R(j) = true;
          continue;
        }
      }
    
      // si on est arrivés jusque là c"est qu"il y a discordance sans remède
      K(i) = false;
      K(j) = false;
      R(j) = false;
      nas++;
    }
  }
  List LL;
  LL["keep"] = K;
  LL["flip"] = R;
  LL["swap_reference"] = ref;
  LL["flip_strand"] = strand;
  LL["NAs"] = nas;
  return LL;
}


RcppExport SEXP gg_alleles_duplicated(SEXP snpsSEXP, SEXP DSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< DataFrame >::type snps(snpsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type D(DSEXP);
    __result = Rcpp::wrap(alleles_duplicated(snps, D));
    return __result;
END_RCPP
}
