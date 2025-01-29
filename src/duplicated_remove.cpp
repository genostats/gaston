#include <Rcpp.h>
#include <iostream>
#include "matrix4.h"
#include "flip_strand.h"

using namespace Rcpp;

inline uint8_t compare_geno(uint8_t x1, uint8_t x2, bool f, bool na, LogicalVector & I, int k) {
  if( f && x2!=3 ) x2=2-x2;
  if( !na && x1==3 ) return x2;
  if( !na && x2==3 ) return x1;
  if(x1!=x2) {
	// Noter incompatibilité
	I(k) = true;
    return 3;
  }
  return(x1);
}

//[[Rcpp::export]]
XPtr<matrix4> duplicated_remove(XPtr<matrix4> x, NumericVector D, LogicalVector keep, LogicalVector flip, int newm, bool na, bool incomp) {
  int n = x->ncol;
  int m = x->nrow;
  
  XPtr<matrix4> r(new matrix4(newm,n));
  int col=0;
  for(int i = 0; i < m; i++) {
    if ( keep(i) == FALSE ) continue;
    if ( keep(i) == TRUE ) {
		
      for(int k = 0; k < n; k++) 
        (*r)(col,k)=(*x)(i,k);
  
      if ( R_IsNA(D(i)) ) {
        col++;
        continue;
      } else {
		// Vecteur pour prendre en compte si il y a eu une incompatibilité avant
		// Afin d'éviter de remettre un génotype si on essaye de remplir les trous
		LogicalVector I(n, false);
		
        for(int j = 0; j < m; j++) {
          if( D(i)!=D(j) || i==j ) {
            continue;
          } else {
            for(int k = 0; k < n; k++)
			  if (!I(k)) (*r)(col,k) = compare_geno( (*r)(col,k), (*x)(j,k), flip(j), na, I, k );
          }
        }
      }
      col++;
    }
	// Incompatibilité
	if ( LogicalVector::is_na(keep(i)) ) {
      if (incomp) continue;
	  if (!incomp) {
		for(int k = 0; k < n; k++) 
		  (*r)(col,k)=(*x)(i,k);

        col++;
	  }
    }
  }
  return r;
}
