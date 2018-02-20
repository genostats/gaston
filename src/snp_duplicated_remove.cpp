#include <Rcpp.h>
#include <iostream>
#include "matrix4.h"
#include "flip_strand.h"

using namespace Rcpp;

inline uint8_t compare_geno(uint8_t x1, uint8_t x2, bool f, bool na) {
  if( f && x2!=3 ) x2=2-x2;
  if( !na & x1==3 ) return x2;
  if( !na & x2==3 ) return x1;
  if(x1!=x2) return 3;
  return(x1);
}

//[[Rcpp::export]]
XPtr<matrix4> duplicated_remove(XPtr<matrix4> x, NumericVector D, LogicalVector keep, LogicalVector flip, int newm, bool na) {
  int n = x->ncol;
  int m = x->nrow;
  NumericVector snp(m);

  XPtr<matrix4> r(new matrix4(newm,n));
  int col=0;
  for(int i = 0; i < m; i++) {
	if ( !keep(i) )
    {
	  continue;
	}
	if (keep(i)) {
	  for(int k = 0; k < n; k++) snp(k)=(*x)(i,k);
	
	  for(int j = 0; j < m; j++) {
	    if( D(i)!=D(j) ) {
		  continue;
		} else {
	      for(int k = 0; k < n; k++) snp(k) = compare_geno( snp(k), (*x)(j,k), flip(j), na );
        }
	  }

      //copie dans la nouvelle bed.matrix
	  for(int k = 0; k < n; k++) (*r)(col,k)=snp(k);
      col++;
	}
  }
  return r;
}

RcppExport SEXP gg_duplicated_remove(SEXP xSEXP, SEXP DSEXP, SEXP keepSEXP, SEXP flipSEXP, SEXP remSEXP, SEXP naSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< XPtr<matrix4> >::type x(xSEXP);
	Rcpp::traits::input_parameter< NumericVector >::type D(DSEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type keep(keepSEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type flip(flipSEXP);
    Rcpp::traits::input_parameter< int >::type rem(remSEXP);
    Rcpp::traits::input_parameter< bool >::type na(naSEXP);
    __result = Rcpp::wrap(duplicated_remove(x, D, keep, flip, rem, na));
    return __result;
END_RCPP
}

