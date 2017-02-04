#include <Rcpp.h>
#include <iostream>
#include "matrix4.h"

using namespace Rcpp;

inline int onetab_str_token(char * & a, char * & token) {
  if(*a == 0) return 0;
  token = a;
  while(*a != 9 && *a != 0) a++; 
  if(*a == 9) {
    *a = 0;
    a++;
    return (a-token-1);
  }  
  return (a-token);
}

// [[Rcpp::export]]
List read_vcf(Function f, int nsamples, int nsnps) {
 
  List L;
  std::vector<std::string> id, ref, alt;
  std::vector<int> chr, pos;

  int chr_, pos_;
  std::string id_("(no SNP read yet)"), ref_, alt_;

  XPtr<matrix4> pX(new matrix4(nsnps, nsamples));
 
  int i = 0;
  while(i < nsnps) {
    SEXP fres = f();
    if (TYPEOF(fres) == LGLSXP && !as<bool>(fres)) {
      break;
    }
    if (TYPEOF(fres) != STRSXP) {
      Rf_error("VCF error, unexpected result from read function");
    }
    CharacterVector x(fres);
    char * a = (char *) x[0];
    char * t = a;

    if(onetab_str_token(a,t)>0) chr_ = atoi(t);
    else Rf_error("VCF format error, last SNP read %s", id_.c_str());

    if(onetab_str_token(a,t)>0) pos_ = atoi(t);
    else Rf_error("VCF format error, last SNP read %s", id_.c_str());

    if(onetab_str_token(a,t)>0) id_.assign(t);
    else Rf_error("VCF format error, last SNP read %s", id_.c_str());

    if(onetab_str_token(a,t)>0) ref_.assign(t);
    else Rf_error("VCF format error while reading SNP read %s", id_.c_str());

    if(onetab_str_token(a,t)>0) alt_.assign(t);
    else Rf_error("VCF format error while reading SNP read %s", id_.c_str());

    if(alt_.find(',') != std::string::npos) continue;
    chr.push_back(chr_);
    pos.push_back(pos_);
    id.push_back(id_);
    ref.push_back(ref_);
    alt.push_back(alt_);

    if(!(onetab_str_token(a,t)>0)) Rf_error("VCF format error while reading SNP read %s", id_.c_str()); // skip qual
    if(!(onetab_str_token(a,t)>0)) Rf_error("VCF format error while reading SNP read %s", id_.c_str()); // skip filter 
    if(!(onetab_str_token(a,t)>0)) Rf_error("VCF format error while reading SNP read %s", id_.c_str()); // skip info
    if(!(onetab_str_token(a,t)>0)) Rf_error("VCF format error while reading SNP read %s", id_.c_str()); // skip format [should check if GT !!]

    for(int j = 0; j < nsamples; j++) {
      int le = onetab_str_token(a,t); 
      if(le != 3 && le != 1) {
        // Rf_error("VCF format error while reading SNP read %s", id_.c_str());
        (*pX).set(i,j,3); // set to NA
      } else if(le == 3) { // deux allèles
        int g = 0;
        if(*t == 49) g++;
        if(*(t+2) == 49) g++;
        if( (*t == 46) || (*(t+2) == 46)) g = 3;
        (*pX).set(i,j,g); 
      } else { // un allèle
        int g = 0;
        if(*t == 49) g = 1; else if(*t == 46) g = 3;
        (*pX).set(i,j,g);
      } 
    }
    i++;
  }

  pX->nrow = i; // si on en a lu moins que prévu !!

  L["id"] = id;
  L["pos"] = pos;
  L["chr"] = chr;
  L["A1"] = ref;
  L["A2"] = alt;
  L["bed"] = pX;
  return L;
}

// [[Rcpp::export]]
int count_dia_vcf(Function f) {
 
  int i = 0;
  std::string alt_;
  char * id_ = const_cast<char *>("(no SNP read yet)");

  for(;;) {
    SEXP fres = f();
    if (TYPEOF(fres) == LGLSXP && !as<bool>(fres)) {
      break;
    }
    if (TYPEOF(fres) != STRSXP) {
      Rf_error("VCF error, unexpected result from read function");
    }
    CharacterVector x(fres);
    char * a = (char *) x[0];
    char * t = a;

    if(!(onetab_str_token(a,t)>0)) Rf_error("VCF format error, last SNP read %s",id_); // skip chr
    if(!(onetab_str_token(a,t)>0)) Rf_error("VCF format error, last SNP read %s",id_); // skip pos

    if(onetab_str_token(a,t)>0) id_ = t ;  // read id
    else Rf_error("VCF format error, last SNP read %s",id_);

    if(!(onetab_str_token(a,t)>0)) Rf_error("VCF format error while reading SNP %s",id_); // skip ref

    if(onetab_str_token(a,t)>0) alt_.assign(t); // read alt
    else Rf_error("VCF format error while reading SNP %s",id_);

    if(alt_.find(',') != std::string::npos) continue;
    i++;
  }

  return i;
}


RcppExport SEXP gg_read_vcf(SEXP fSEXP, SEXP nsamplesSEXP, SEXP nsnpsSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< Function >::type f(fSEXP );
        Rcpp::traits::input_parameter< int >::type nsamples(nsamplesSEXP );
        Rcpp::traits::input_parameter< int >::type nsnps(nsnpsSEXP );
        List __result = read_vcf(f, nsamples, nsnps);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}


RcppExport SEXP gg_count_dia_vcf(SEXP fSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< Function >::type f(fSEXP );
        int __result = count_dia_vcf(f);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}

