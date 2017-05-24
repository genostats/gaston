#include <Rcpp.h>
#include <iostream>
#include "matrix4.h"
#include <fstream>
#include "gzstream.h"

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


// skip header and get samples
CharacterVector vcf_samples(igzstream & in) {
  std::string line;
  std::vector<std::string> sa;
  while(std::getline(in, line)) {
    if(line.substr(0,1) != "#") stop("Bad VCF format");
    if(line.substr(0,2) != "##") break; // fin des meta information
  }
  char * a = &line[0];
  char * t = a;
  for(int i = 0; i < 9; i++) {
    if(str_token_tab(a, t) <= 0) stop("Bad VCF format");
  }
  while(str_token_tab(a, t) > 0)
    sa.push_back( std::string(t) );
  return wrap(sa);
}



List read_vcf2(CharacterVector filename, int max_snps, bool get_info) {

  if(filename.length() != 1) stop("filename should be a CharacterVector of length 1");
  igzstream in( (char *) filename[0] );
  if(!in.good()) stop("Can't open file");
  CharacterVector samples = vcf_samples(in);
  int nsamples = samples.length();

  List L;
  std::vector<std::string> chr, id, ref, alt, filter, info;
  std::vector<int> pos;
  std::vector<double> qual;

  XPtr<matrix4> pX(new matrix4(0, nsamples));  // avec nrow = 0 allocations() n'est pas appelé

  std::vector<uint8_t *> data;
  uint8_t * data_;

  std::string line;
  int i = 0;
  while(i != max_snps && std::getline(in, line)) {
    int pos_ = 0;
    std::string chr_, id_("(no SNP read yet)"), ref_, alt_, filter_, info_;
    double qual_;

    char * a = &line[0];
    char * t = a;

    if(str_token_tab(a,t)>0) chr_.assign(t);     // CHR
    else Rf_error("VCF format error, last SNP read %s chr = %s pos %d", id_.c_str(), chr_.c_str(), pos_);

    if(str_token_tab(a,t)>0) pos_ = atoi(t);     // POS
    else Rf_error("VCF format error, last SNP read %s chr = %s pos %d", id_.c_str(), chr_.c_str(), pos_);

    if(str_token_tab(a,t)>0) id_.assign(t);      // ID
    else Rf_error("VCF format error, last SNP read %s chr = %s pos %d", id_.c_str(), chr_.c_str(), pos_);

    if(str_token_tab(a,t)>0) ref_.assign(t);     // REF
    else Rf_error("VCF format error while reading SNP %s chr = %s pos %d", id_.c_str(), chr_.c_str(), pos_);

    if(str_token_tab(a,t)>0) alt_.assign(t);     // ALT
    else Rf_error("VCF format error while reading SNP %s chr = %s pos %d", id_.c_str(), chr_.c_str(), pos_);

    if(alt_.find(',') != std::string::npos) continue;  // on ne continue que s'il y a un seul allèle alternatif

    if(str_token_tab(a,t)>0) qual_ = atof(t);    // QUAL
    else Rf_error("VCF format error while reading SNP %s chr = %s pos %d", id_.c_str(), chr_.c_str(), pos_);

    if(str_token_tab(a,t)>0) filter_.assign(t);  // FILTER
    else Rf_error("VCF format error while reading SNP %s chr = %s pos %d", id_.c_str(), chr_.c_str(), pos_);

    if(str_token_tab(a,t)>0) info_.assign(t);    // INFO
    else Rf_error("VCF format error while reading SNP %s chr = %s pos %d", id_.c_str(), chr_.c_str(), pos_);

    if(str_token_tab(a,t)>0) {                   // FORMAT
      char * b;
      if(str_token_col(t,b)>0) {
        if(strcmp(b,"GT") != 0) continue;        // doit commencer par GT.
      }
      else Rf_error("VCF format error while reading SNP %s chr = %s pos %d", id_.c_str(), chr_.c_str(), pos_);
    } else Rf_error("VCF format error while reading SNP %s chr = %s pos %d", id_.c_str(), chr_.c_str(), pos_);


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

    int j = 0;
    for(int j1 = 0; j1 < nsamples; j1++) {
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
      j++;
    }

    data.push_back(data_);
    i++;
  }
  
  // et on finit la construction de la matrice
  pX->nrow = i; 
  if(i > 0) {
    pX->data = new uint8_t * [i];
    for(int j = 0; j < i; j++) pX->data[j] = data[j];
  }

  L["id"] = id;
  L["pos"] = pos;
  L["chr"] = chr;
  L["A1"] = ref;
  L["A2"] = alt;
  L["quality"] = qual;
  L["filter"] = filter;
  if(get_info) L["info"] = info;

  L["bed"] = pX;
  L["samples"] = samples;
  return L;
}


// lecture VCF filtré par chr / pos
// POS est une liste de positions dont les clefs sont les ids des chromosomes
// on suppose que les positions sont triées
// aucune supposition n'est faite sur l'ordre des variants dans le VCF
inline bool is_ok(std::string chr, int pos, List POS) {
  if(!POS.containsElementNamed(chr.c_str())) {
    return false;
  }
  IntegerVector a( as<IntegerVector>(POS[chr.c_str()]) );
  return std::binary_search(a.begin(), a.end(), pos);
}

// [[Rcpp::export]]
List read_vcf_filtered(CharacterVector filename, List POS, int max_snps, bool get_info) {
 
  if(filename.length() != 1) stop("filename should be a CharacterVector of length 1");
  igzstream in( (char *) filename[0] );
  if(!in.good()) stop("Can't open file");
  CharacterVector samples = vcf_samples(in);
  int nsamples = samples.length();

  List L;
  std::vector<std::string> chr, id, ref, alt, filter, info;
  std::vector<int> pos;
  std::vector<double> qual;
 
  XPtr<matrix4> pX(new matrix4(0, nsamples));  // avec nrow = 0 allocations() n'est pas appelé

  std::vector<uint8_t *> data;
  uint8_t * data_;

  std::string line;
  int i = 0;
  while(i != max_snps && std::getline(in, line)) {
    int pos_ = 0;
    std::string chr_, id_("(no SNP read yet)"), ref_, alt_, filter_, info_;
    double qual_;

    char * a = &line[0];
    char * t = a;

    if(str_token_tab(a,t)>0) chr_.assign(t);     // CHR
    else Rf_error("VCF format error, last SNP read %s chr = %s pos %d", id_.c_str(), chr_.c_str(), pos_);

    if(str_token_tab(a,t)>0) pos_ = atoi(t);     // POS
    else Rf_error("VCF format error, last SNP read %s chr = %s pos %d", id_.c_str(), chr_.c_str(), pos_);

    if(!is_ok(chr_, pos_, POS)) continue;        // filtrage

    if(str_token_tab(a,t)>0) id_.assign(t);      // ID
    else Rf_error("VCF format error, last SNP read %s chr = %s pos %d", id_.c_str(), chr_.c_str(), pos_);

    if(str_token_tab(a,t)>0) ref_.assign(t);     // REF
    else Rf_error("VCF format error while reading SNP %s chr = %s pos %d", id_.c_str(), chr_.c_str(), pos_);

    if(str_token_tab(a,t)>0) alt_.assign(t);     // ALT
    else Rf_error("VCF format error while reading SNP %s chr = %s pos %d", id_.c_str(), chr_.c_str(), pos_);

    if(alt_.find(',') != std::string::npos) continue;  // on ne continue que s'il y a un seul allèle alternatif

    if(str_token_tab(a,t)>0) qual_ = atof(t);    // QUAL
    else Rf_error("VCF format error while reading SNP %s chr = %s pos %d", id_.c_str(), chr_.c_str(), pos_);

    if(str_token_tab(a,t)>0) filter_.assign(t);  // FILTER
    else Rf_error("VCF format error while reading SNP %s chr = %s pos %d", id_.c_str(), chr_.c_str(), pos_);

    if(str_token_tab(a,t)>0) info_.assign(t);    // INFO
    else Rf_error("VCF format error while reading SNP %s chr = %s pos %d", id_.c_str(), chr_.c_str(), pos_);


    if(str_token_tab(a,t)>0) {                   // FORMAT
      char * b;
      if(str_token_col(t,b)>0) {
        if(strcmp(b,"GT") != 0) continue;        // doit commencer par GT.
      }
      else Rf_error("VCF format error while reading SNP %s chr = %s pos %d", id_.c_str(), chr_.c_str(), pos_);
    } else Rf_error("VCF format error while reading SNP %s chr = %s pos %d", id_.c_str(), chr_.c_str(), pos_);


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

    int j = 0;
    for(int j1 = 0; j1 < nsamples; j1++) {
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
      j++;
    }

    data.push_back(data_);
    i++;
  }
  
  // et on finit la construction de la matrice
  pX->nrow = i; 
  if(i > 0) {
    pX->data = new uint8_t * [i];
    for(int j = 0; j < i; j++) pX->data[j] = data[j];
  }

  L["id"] = id;
  L["pos"] = pos;
  L["chr"] = chr;
  L["A1"] = ref;
  L["A2"] = alt;
  L["quality"] = qual;
  L["filter"] = filter;
  if(get_info) L["info"] = info;

  L["bed"] = pX;
  L["samples"] = samples;
  return L;
}

RcppExport SEXP gg_read_vcf2(SEXP fSEXP, SEXP maxsnpSEXP, SEXP giSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< CharacterVector >::type f(fSEXP );
        Rcpp::traits::input_parameter< int >::type maxsnp(maxsnpSEXP );
        Rcpp::traits::input_parameter< bool >::type gi(giSEXP );
        List __result = read_vcf2(f, maxsnp, gi);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}

RcppExport SEXP gg_read_vcf_filtered(SEXP fSEXP, SEXP POSSEXP, SEXP max_snpsSEXP, SEXP get_infoSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type f(fSEXP);
    Rcpp::traits::input_parameter< List >::type POS(POSSEXP);
    Rcpp::traits::input_parameter< int >::type max_snps(max_snpsSEXP);
    Rcpp::traits::input_parameter< bool >::type get_info(get_infoSEXP);
    rcpp_result_gen = Rcpp::wrap(read_vcf_filtered(f, POS, max_snps, get_info));
    return rcpp_result_gen;
END_RCPP
}

