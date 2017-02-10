#include <Rcpp.h>
#include <iostream>
#include "matrix4.h"

using namespace Rcpp;

// renvoie la longueur du token trouvé
// token = a après l'appel ; on peut le lire directement
// a a été modifié pour être placé après le délimiteur trouvé qu'on a remplacé par 0
// du coup on peut rappeler la même fonction
inline int str_token(char * & a, char * & token, char delim) {
  if(*a == 0) return 0;
  token = a;
  while(*a != delim && *a != 0) a++; 
  if(*a == delim) {
    *a = 0;
    a++;
    return (a-token-1);
  }  
  return (a-token);
}

inline int str_token_tab(char * & a, char * & token) {
  return str_token(a, token, 9);
}

void set(uint8_t * data, size_t j, uint8_t val) {
  uint8_t & a = data[j/4];
  a &= ~(3 << ((j%4)*2));  // set to 00
  a |= (val << ((j%4)*2)); // set to val
}


// [[Rcpp::export]]
List read_vcf2(Function f, int nsamples, bool get_info) {
 
  List L;
  std::vector<std::string> id, ref, alt, filter, info;
  std::vector<int> chr, pos;
  std::vector<double> qual;

  int chr_, pos_;
  std::string id_("(no SNP read yet)"), ref_, alt_, filter_, info_;

  XPtr<matrix4> pX(new matrix4(0, nsamples));  // avec nrow = 0 allocations() n'est pas appelé

  std::vector<uint8_t *> data;
  uint8_t * data_;

  int i = 0;
  while(true) {
    SEXP fres = f();
    if (TYPEOF(fres) == LGLSXP && !as<bool>(fres)) { // eof
      break;
    }
    if (TYPEOF(fres) != STRSXP) {
      Rf_error("VCF error, unexpected result from read function");
    }
    CharacterVector x(fres);
    char * a = (char *) x[0];
    char * t = a;

    if(str_token_tab(a,t)>0) chr_ = atoi(t);     // CHR
    else Rf_error("VCF format error, last SNP read %s", id_.c_str());

    if(str_token_tab(a,t)>0) pos_ = atoi(t);     // POS
    else Rf_error("VCF format error, last SNP read %s", id_.c_str());

    if(str_token_tab(a,t)>0) id_.assign(t);      // ID
    else Rf_error("VCF format error, last SNP read %s", id_.c_str());

    if(str_token_tab(a,t)>0) ref_.assign(t);     // REF
    else Rf_error("VCF format error while reading SNP read %s", id_.c_str());

    if(str_token_tab(a,t)>0) alt_.assign(t);     // ALT
    else Rf_error("VCF format error while reading SNP read %s", id_.c_str());

    if(alt_.find(',') != std::string::npos) continue;  // un seul allèle alternatif

    // maintenant qu'on sait qu'on va garder cette ligne, on fait nos push back
    chr.push_back(chr_);
    pos.push_back(pos_);
    id.push_back(id_);
    ref.push_back(ref_);
    alt.push_back(alt_);

    if(str_token_tab(a,t)>0) qual.push_back(atof(t));  // QUAL
    else Rf_error("VCF format error while reading SNP read %s", id_.c_str()); 

    if(str_token_tab(a,t)>0) {       // FILTER
      filter_.assign(t);
      filter.push_back(filter_);
    } else Rf_error("VCF format error while reading SNP read %s", id_.c_str());

    if(str_token_tab(a,t)>0) {      // INFO
      if(get_info) {
        info_.assign(t);
        info.push_back(info_);
      }
    } else Rf_error("VCF format error while reading SNP read %s", id_.c_str()); 

    if(!(str_token_tab(a,t)>0)) Rf_error("VCF format error while reading SNP read %s", id_.c_str()); // skip format [should check if GT !!]

    // nouvelle ligne de données    
    data_ = new uint8_t [pX->true_ncol];
    std::fill(data_, data_ + pX->true_ncol, 255); // c'est important de remplir avec 3 -> NA

    for(int j = 0; j < nsamples; j++) {
      int le = str_token_tab(a,t); 
      if(le != 3 && le != 1) {
        // Rf_error("VCF format error while reading SNP read %s", id_.c_str());
        // (*pX).set(i,j,3); // set to NA
        set(data_, j, 3);
      } else if(le == 3) { // deux allèles
        int g = 0;
        if(*t == 49) g++;
        if(*(t+2) == 49) g++;
        if( (*t == 46) || (*(t+2) == 46)) g = 3;
        // (*pX).set(i,j,g); 
        set(data_, j, g);
      } else { // un allèle
        int g = 0;
        if(*t == 49) g = 1; else if(*t == 46) g = 3;
        // (*pX).set(i,j,g);
        set(data_, j, g);
      } 
    }

    data.push_back(data_);
    i++;
  }
  
  // et on finit la construction de la matrice
  pX->nrow = i; 
  pX->data = new uint8_t * [i];
  for(int j = 0; j < i; j++) pX->data[j] = data[j];

  L["id"] = id;
  L["pos"] = pos;
  L["chr"] = chr;
  L["A1"] = ref;
  L["A2"] = alt;
  L["quality"] = qual;
  L["filter"] = filter;
  if(get_info) L["info"] = info;

  L["bed"] = pX;
  return L;
}

RcppExport SEXP gg_read_vcf2(SEXP fSEXP, SEXP nsamplesSEXP, SEXP giSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< Function >::type f(fSEXP );
        Rcpp::traits::input_parameter< int >::type nsamples(nsamplesSEXP );
        Rcpp::traits::input_parameter< bool >::type gi(giSEXP );
        List __result = read_vcf2(f, nsamples, gi);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}

