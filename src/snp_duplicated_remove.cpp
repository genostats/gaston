#include <Rcpp.h>
#include <iostream>
#include "matrix4.h"
#include "flip_strand.h"

using namespace Rcpp;

inline uint8_t compare_geno(uint8_t x1, uint8_t x2, std::string f) {
  if(f=="ref" && x2!=3) x2=2-x2;
  if(x1!=x2) return 3;
  return(x1);
}

//[[Rcpp::export]]
XPtr<matrix4> duplicated_remove(XPtr<matrix4> x, NumericVector D, StringVector flip, int rem) {
  int n = x->ncol;
  int m = x->nrow;
  NumericVector snp(m);
  int nn=n;
  nn -= rem;

  XPtr<matrix4> r(new matrix4(m,nn));
  int col=0;
  for(int i = 0; i < n; i++) {
	if ( flip(i)=="remove" || flip(i)=="ref")
    {
	  continue;
	}
	if (flip(i)=="keep") {
	  for(int k = 0; k < m; k++) snp(k)=(*x)(i,k);
	
	  for(int j = 0; j < n; j++) {
	    if( D(i)!=D(j) ) {
		  continue;
		} else {
	      for(int k = 0; k < m; k++) snp(k) = compare_geno( snp(k), (*x)(j,k), Rcpp::as< std::string >( flip(j) ) );
        }
	  }

      //copie dans la nouvelle bed.matrix
	  for(int k = 0; k < m; k++) (*r)(col++,k)=snp(k);
	}
  }
  return r;
}

RcppExport SEXP gg_duplicated_remove(SEXP xSEXP, SEXP DSEXP, SEXP flipSEXP, SEXP remSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< XPtr<matrix4> >::type x(xSEXP);
	Rcpp::traits::input_parameter< NumericVector >::type D(DSEXP);
    Rcpp::traits::input_parameter< StringVector >::type flip(flipSEXP);
    Rcpp::traits::input_parameter< int >::type rem(remSEXP);
    __result = Rcpp::wrap(duplicated_remove(x, D, flip, rem));
    return __result;
END_RCPP
}

