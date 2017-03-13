#include <Rcpp.h>
#include "matrix-varia.h"
#include "matrix4.h"

List GWAS_lmm_score(XPtr<matrix4> pA, NumericVector PY, NumericMatrix P, NumericVector mu, int beg, int end) {
  Map_MatrixXd Py(as<Map<MatrixXd> >(PY));
  Map_MatrixXd PP(as<Map<MatrixXd> >(P));

  int r= end-beg+1;
  int n=Py.rows();
  VectorXd SNP(n);
  NumericVector s(r);
  double t, v;  
  
  for(int i = beg; i <= end; i++) {
    if( std::isnan(mu(i)) || mu(i) == 0 || mu(i) == 2 ) {
      s(i-beg) = NAN;
      continue;
    }
    // récupérer SNP
    for(int ii = 0; ii < pA->true_ncol-1; ii++) {
      uint8_t x = pA->data[i][ii];
      for(int ss = 0; ss < 4; ss++) {
        SNP(4*ii+ss) = ((x&3) != 3)?(x&3):mu(i);
        x >>= 2;
      }
    }
    { int ii = pA->true_ncol-1;
      uint8_t x = pA->data[i][ii];
      for(int ss = 0; ss < 4 && 4*ii+ss < pA->ncol; ss++) {
        SNP(4*ii+ss) = ((x&3) != 3)?(x&3):mu(i);
        x >>= 2;
      }
    }
    
    v = (PP.selfadjointView<Lower>()*SNP).dot(SNP);
    t = SNP.dot(Py.col(0));
    s(i-beg) = t*t/v;
  }
  
  List S;
  S["score"] = s;

  return S;
}


// float version
//[[Rcpp::export]]
List GWAS_lmm_score_f(XPtr<matrix4> pA, NumericVector PY, NumericMatrix P, NumericVector mu, int beg, int end) {
  int n = PY.size();
  if(P.nrow() != n || P.ncol() != n) 
    stop("Dimensions mismatch\n");

  VectorXf Py(n);
  for(int i = 0; i < n; i++) Py(i) = PY[i];

  MatrixXf PP(n,n);
  for(int i = 0; i < n; i++)
    for(int j = 0; j < n; j++) 
      PP(i,j) = P(i,j);
     

  int r = end-beg+1;

  VectorXf SNP(n);
  NumericVector s(r);
  double t, v;  
  
  for(int i = beg; i <= end; i++) {
    if( std::isnan(mu(i)) || mu(i) == 0 || mu(i) == 2 ) {
      s(i-beg) = NAN;
      continue;
    }
    // récupérer SNP
    for(int ii = 0; ii < pA->true_ncol-1; ii++) {
      uint8_t x = pA->data[i][ii];
      for(int ss = 0; ss < 4; ss++) {
        SNP(4*ii+ss) = ((x&3) != 3)?(x&3):mu(i);
        x >>= 2;
      }
    }
    { int ii = pA->true_ncol-1;
      uint8_t x = pA->data[i][ii];
      for(int ss = 0; ss < 4 && 4*ii+ss < pA->ncol; ss++) {
        SNP(4*ii+ss) = ((x&3) != 3)?(x&3):mu(i);
        x >>= 2;
      }
    }
    
    v = (PP.selfadjointView<Lower>()*SNP).dot(SNP);
    t = SNP.dot(Py);
    s(i-beg) = t*t/v;
  }
  
  List S;
  S["score"] = s;

  return S;
}


RcppExport SEXP gg_GWAS_lmm_score(SEXP pASEXP, SEXP PYSEXP, SEXP PSEXP, SEXP muSEXP, SEXP begSEXP, SEXP endSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< XPtr<matrix4> >::type pA(pASEXP );
        Rcpp::traits::input_parameter< NumericVector >::type PY(PYSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type P(PSEXP );
	Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP );
        Rcpp::traits::input_parameter< int >::type beg(begSEXP );
        Rcpp::traits::input_parameter< int >::type end(endSEXP );
        List __result = GWAS_lmm_score(pA, PY, P, mu, beg, end);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}

RcppExport SEXP gg_GWAS_lmm_score_f(SEXP pASEXP, SEXP PYSEXP, SEXP PSEXP, SEXP muSEXP, SEXP begSEXP, SEXP endSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< XPtr<matrix4> >::type pA(pASEXP);
    Rcpp::traits::input_parameter< NumericVector >::type PY(PYSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type P(PSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< int >::type beg(begSEXP);
    Rcpp::traits::input_parameter< int >::type end(endSEXP);
    rcpp_result_gen = Rcpp::wrap(GWAS_lmm_score_f(pA, PY, P, mu, beg, end));
    return rcpp_result_gen;
END_RCPP
}

