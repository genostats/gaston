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

inline int str_token_col(char * & a, char * & token) {  
  return str_token(a, token, ':');
}

void set(uint8_t * data, size_t j, uint8_t val) {
  uint8_t & a = data[j/4];
  a &= ~(3 << ((j%4)*2));  // set to 00
  a |= (val << ((j%4)*2)); // set to val
}


// [[Rcpp::export]]
List read_vcf2(Function f, int nsamples, int max_snps, bool get_info) {
 
  List L;
  std::vector<std::string> id, ref, alt, filter, info;
  std::vector<int> chr, pos;
  std::vector<double> qual;
  
  XPtr<matrix4> pX(new matrix4(0, nsamples));  // avec nrow = 0 allocations() n'est pas appelé

  std::vector<uint8_t *> data;
  uint8_t * data_;

  int i = 0;
  while(i != max_snps) {
    int chr_, pos_;
    std::string id_("(no SNP read yet)"), ref_, alt_, filter_, info_;
    double qual_;

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

    if(alt_.find(',') != std::string::npos) continue;  // on ne continue que s'il y a un seul allèle alternatif

    if(str_token_tab(a,t)>0) qual_ = atof(t);    // QUAL
    else Rf_error("VCF format error while reading SNP read %s", id_.c_str()); 

    if(str_token_tab(a,t)>0) filter_.assign(t);  // FILTER
    else Rf_error("VCF format error while reading SNP read %s", id_.c_str());

    if(str_token_tab(a,t)>0) info_.assign(t);    // INFO
    else Rf_error("VCF format error while reading SNP read %s", id_.c_str()); 

    if(str_token_tab(a,t)>0) {                   // FORMAT
      char * b;
      if(str_token_col(t,b)>0) {
        if(strcmp(b,"GT") != 0) continue;        // doit commencer par GT.
      }
      else Rf_error("VCF format error while reading SNP read %s", id_.c_str());
    } else Rf_error("VCF format error while reading SNP read %s", id_.c_str());


    // maintenant qu'on sait qu'on va garder cette ligne, on fait nos push back
    chr.push_back(chr_);
    pos.push_back(pos_);
    id.push_back(id_);
    ref.push_back(ref_);
    alt.push_back(alt_);
    qual.push_back(qual_); 
    filter.push_back(filter_);
    if(get_info) info.push_back(info_);
    // nouvelle ligne de données    
    data_ = new uint8_t [pX->true_ncol];
    std::fill(data_, data_ + pX->true_ncol, 255); // c'est important de remplir avec 3 -> NA

    for(int j = 0; j < nsamples; j++) {
      char * b;
      int g = 0;
      if(str_token_tab(a,t) == 0)
        Rf_error("VCF format error while reading SNP read %s", id_.c_str());
      int le = str_token_col(t,b);
      if(le == 3) { // deux allèles 0/0 0/1 1/1 ou 0|0 etc  
        if(*b == '1') g++;
        if(*(b+2) == '1') g++;
        if( (*b == '.') || (*(b+2) == '.')) g = 3;  // valeurs manquantes = .
      } else if(le == 1) { // un allèle (cas haploide)
        if(*b == '1') g = 1; else if(*b == '.') g = 3;
      } else {
        g = 3; // set to NA
      }
      set(data_, j, g);
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

RcppExport SEXP gg_read_vcf2(SEXP fSEXP, SEXP nsamplesSEXP, SEXP maxsnpSEXP, SEXP giSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< Function >::type f(fSEXP );
        Rcpp::traits::input_parameter< int >::type nsamples(nsamplesSEXP );
        Rcpp::traits::input_parameter< int >::type maxsnp(maxsnpSEXP );
        Rcpp::traits::input_parameter< bool >::type gi(giSEXP );
        List __result = read_vcf2(f, nsamples, maxsnp, gi);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}

