#include <Rcpp.h>
#include <iostream>
#include "flip_strand.h"

using namespace Rcpp;

// ci-dessous il faut vérifier en amont que
// * tous les data frames ont des composantes A1 et A2
// * il n"y a que des allèles A C G T

// [[Rcpp::export]]
List alleles_duplicated(DataFrame snps, NumericVector D) {
  int n = snps.nrows();

  int nas = 0, strand = 0, ref = 0;
  
  StringVector R(n);
  CharacterVector A1 = as<CharacterVector>(snps["A1"]);
  CharacterVector A2 = as<CharacterVector>(snps["A2"]);
  
  for(int i = 0; i < n; i++) { // on parcourt les SNP
  
    if ( R_IsNA(D(i)) || R(i)=="ref" || R(i)=="remove" ) continue; // pas de dupliqué ou déjà traité
	
	R(i) = "keep";
	
    const char * a1 = CHAR(STRING_ELT(A1,i));
    const char * a2 = CHAR(STRING_ELT(A2,i));
	
	for(int j=i+1; j<n; j++) { // on cherche les duplica
      if( D(i)!=D(j) ) {
        continue;
	  } else {
        const char * b1 = CHAR(STRING_ELT(A1,j));
        const char * b2 = CHAR(STRING_ELT(A2,j));
        // tout ok
        if( !std::strcmp(a1, b1) && !std::strcmp(a2, b2) ) {
          R(j) = "remove";
          continue; // case closed
        }
        // inversion allèle ref ?
        if( !std::strcmp(a1, b2) && !std::strcmp(a2, b1) ) {
          R(j) = "ref";
          ref++;
          continue; // case closed
        }
        // inversion de brin
        std::string c1 = flip_strand(b1);
        std::string c2 = flip_strand(b2);
        if( c1 == a1 && c2 == a2 ) { // ok
          R(j) = "remove";
          strand++;
          continue;
        }
        if( c1 == a2 && c2 == a1 ) { // + inversion allèle ref
          R(j) = "ref";
          strand++; 
          ref++;
          continue;
        }
		// un des cas monomorphes
		if( !std::strcmp(a1, "0") ) {
		  if( !std::strcmp(a2, b2) ) {
		    R(i) = "remove";
			R(j) = "keep";
			// Il vaut mieux garder le second qui repassera dans la boucle
			// pour un potentiel autre duplica
			break;
		  }		  
		  if( !std::strcmp(a2, b1) ) {
		    R(i) = "ref";
			R(j) = "keep";
			break;
		  }
		  std::string c2 = flip_strand(b2);
          if( a2==c2 ) {
		     R(i) = "remove";
			 R(j) = "keep";
			 break;
		  }
		  std::string c1 = flip_strand(b1);
          if( a2==c1 ) {
		     R(i) = "ref";
			 R(j) = "keep";
			 break;
		  }
	    }
		if( !std::strcmp(a2, "0") ) {
		  if( !std::strcmp(a1, b1) ) {
		    R(i) = "remove";
			R(j) = "keep";
			break;
		  }		  
		  if( !std::strcmp(a1, b2) ) {
		    R(i) = "ref";
			R(j) = "keep";
			break;
		  }
		  std::string c1 = flip_strand(b1);
          if( a1==c1 ) {
		     R(i) = "remove";
			 R(j) = "keep";
			 break;
		  }
		  std::string c2 = flip_strand(b2);
          if( a1==c2 ) {
		     R(i) = "ref";
			 R(j) = "keep";
			 break;
		  }
	    }
		if( !std::strcmp(b1, "0") ) {
		  if( !std::strcmp(a2, b2) ) {
		    R(i) = "keep";
			R(j) = "remove";
			continue;
		  }		  
		  if( !std::strcmp(a1, b2) ) {
		    R(i) = "keep";
			R(j) = "ref";
			continue;
		  }
		  std::string c2 = flip_strand(a2);
          if( c2==b2 ) {
		     R(i) = "keep";
			 R(j) = "remove";
			 continue;
		  }
		  std::string c1 = flip_strand(a1);
          if( c1==b2 ) {
		     R(i) = "keep";
			 R(j) = "ref";
			 continue;
		  }
	    }
		if( !std::strcmp(b2, "0") ) {
		  if( !std::strcmp(a1, b1) ) {
		    R(i) = "keep";
			R(j) = "remove";
			continue;
		  }		  
		  if( !std::strcmp(a2, b1) ) {
		    R(i) = "keep";
			R(j) = "ref";
			continue;
		  }
		  std::string c1 = flip_strand(a1);
          if( c1==b1 ) {
		     R(i) = "keep";
			 R(j) = "remove";
			 continue;
		  }
		  std::string c2 = flip_strand(a2);
          if( c2==b1 ) {
		     R(i) = "keep";
			 R(j) = "ref";
			 continue;
		  }
	    }
		
        // si on est arrivés jusque là c"est qu"il y a discordance sans remède
	    R(i) = "remove";
		R(j) = "remove";
        nas++;
      }
	}
  }
  List LL;
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
