//#include <iostream>
#include <string>
#include "string_equal.h"
#include "flip_strand.h"

#ifndef GASTON_SNP_HASH
#define GASTON_SNP_HASH
//  inspiré par le code de Dirk Eddelbuettel, Romain Francois and Kevin Ushey
//  pour la fonction match de Rccp sugar
using namespace Rcpp;

enum hash_type { snpid, snpid_chr_pos, snpid_chr_pos_al, chr_pos, chr_pos_al };

class SNPhash {
  public:
    int n, m, k;
    CharacterVector id;
    IntegerVector chr, pos;
    CharacterVector A1, A2;
    hash_type htype;
    std::vector<int> index;
    int nb_duplicates;
    std::vector<int> dup_indices;


  // empty constructor
  SNPhash() : m(0), k(-1) {}

  // creates Hash by chr:pos
  SNPhash(IntegerVector CHR, IntegerVector POS) : m(2), k(1), chr(CHR), pos(POS), htype(chr_pos) {
    n = chr.length();
    if(pos.length() != n) 
      stop("Length mismatch");
    while( m < 2*n ){ m *= 2; k++; }
    index.resize(m);
    std::fill(index.begin(), index.end(), 0);
    nb_duplicates = 0;
    for(int i = 0; i < n; i++) {
      int p = pos[i], c = chr[i];
      unsigned int ad = hash(32*p + c);
      while(index[ad] && (pos[index[ad] - 1] != p || chr[index[ad] -1] != c)) {
        ++ad %= m;
      }
      if(!index[ad]) {
        index[ad] = i+1;
      } else { 
        nb_duplicates++;
        dup_indices.push_back(i+1);
      }
    }
  }

  // creates Hash by chr:pos:alleles
  SNPhash(IntegerVector CHR, IntegerVector POS, CharacterVector AL1, CharacterVector AL2) 
      : m(2), k(1), chr(CHR), pos(POS), A1(AL1), A2(AL2), htype(chr_pos_al) {
    n = chr.length();
    if(pos.length() != n || A1.length() != n || A2.length() != n) 
      stop("Length mismatch");
    while( m < 2*n ){ m *= 2; k++; }
    index.resize(m);
    std::fill(index.begin(), index.end(), 0);
    nb_duplicates = 0;
    for(int i = 0; i < n; i++) {
      int p = pos[i], c = chr[i];
      const char * a1 = CHAR(STRING_ELT(A1,i));
      const char * a2 = CHAR(STRING_ELT(A2,i));
      unsigned int ad = hash(32*p + c);
      while(index[ad] && (pos[index[ad] - 1] != p || chr[index[ad] -1] != c 
                          || std::strcmp(a1, CHAR(STRING_ELT(A1,index[ad] -1))) 
                          || std::strcmp(a2, CHAR(STRING_ELT(A2,index[ad] -1))))) {
        ++ad %= m;
      }
      if(!index[ad]) {
        index[ad] = i+1;
      } else { 
        nb_duplicates++;
        dup_indices.push_back(i+1);
      }
    }
  }

  // creates Hash by id
  SNPhash(CharacterVector ID) : m(2), k(1), id(ID), htype(snpid) {
    n = id.length();
    while( m < 2*n ){ m *= 2; k++; }
    index.resize(m);
    std::fill(index.begin(), index.end(), 0);
    nb_duplicates = 0;
    for(int i = 0; i < n; i++) {
      const char * id_ = CHAR(STRING_ELT(id,i));
      unsigned int ad = djb2(id_);
      while(index[ad] && std::strcmp(id_, CHAR(STRING_ELT(id, index[ad] - 1)))) {
        ++ad %= m;
      }
      if(!index[ad]) {
        index[ad] = i+1;
      } else { 
        nb_duplicates++;
        dup_indices.push_back(i+1);
      }
    }
  }

  // creates Hash by id:chr:pos
  SNPhash(CharacterVector ID, IntegerVector CHR, IntegerVector POS) : m(2), k(1), id(ID), chr(CHR), pos(POS), htype(snpid_chr_pos) {
    n = id.length();
    if(chr.length() != n || pos.length() != n) 
      stop("Length mismatch");
    while( m < 2*n ){ m *= 2; k++; }
    index.resize(m);
    std::fill(index.begin(), index.end(), 0);
    nb_duplicates = 0;
    for(int i = 0; i < n; i++) {
      int p = pos[i], c = chr[i];
      const char * id_ = CHAR(STRING_ELT(id,i));
      unsigned int ad = hash_combine( djb2(id_), hash(32*p + c) );
      while(index[ad] && 
            (pos[index[ad] - 1] != p || chr[index[ad] -1] != c 
              || std::strcmp(id_, CHAR(STRING_ELT(id, index[ad] - 1))))) {
        ++ad %= m;
      }
      if(!index[ad]) {
        index[ad] = i+1;
      } else { 
        nb_duplicates++;
        dup_indices.push_back(i+1);
      }
    }
  }

  // creates Hash by id:chr:pos:alleles
  SNPhash(CharacterVector ID, IntegerVector CHR, IntegerVector POS, CharacterVector AL1, CharacterVector AL2) 
      : m(2), k(1), id(ID), chr(CHR), pos(POS), A1(AL1), A2(AL2), htype(snpid_chr_pos_al) {
    n = id.length();
    if(chr.length() != n || pos.length() != n || A1.length() != n || A2.length() != n) stop("Length mismatch");
    while( m < 2*n ){ m *= 2; k++; }
    index.resize(m);
    std::fill(index.begin(), index.end(), 0);
    nb_duplicates = 0;
    for(int i = 0; i < n; i++) {
      int p = pos[i], c = chr[i];
      const char * id_ = CHAR(STRING_ELT(id,i));
      const char * a1 = CHAR(STRING_ELT(A1,i));
      const char * a2 = CHAR(STRING_ELT(A2,i));
      unsigned int ad = hash_combine( djb2(id_), hash(32*p + c) );
      while(index[ad] && (pos[index[ad] - 1] != p || chr[index[ad] -1] != c 
                          || std::strcmp(id_, CHAR(STRING_ELT(id, index[ad] - 1)))
                          || std::strcmp(a1, CHAR(STRING_ELT(A1,index[ad] -1))) 
                          || std::strcmp(a2, CHAR(STRING_ELT(A2,index[ad] -1))))) {
        ++ad %= m;
      }
      if(!index[ad]) {
        index[ad] = i+1;
      } else { 
        nb_duplicates++;
        dup_indices.push_back(i+1);
      }
    }
  }


  // look up by snpid 
  template<typename T>
  inline unsigned int lookup(T ID) const {
    if(htype != snpid)
      return NA_INTEGER;
    unsigned int ad = djb2(ID);
    while(index[ad]) {
      if(string_equal(ID, id[index[ad] - 1]))
        return index[ad];
      ++ad %= m;
    }
    return NA_INTEGER;
  }

  // look up by snpid:chr:pos
  template<typename T>
  inline unsigned int lookup(T ID, int c, int p) const {
    if(htype != snpid_chr_pos && htype != snpid_chr_pos_al)
      return NA_INTEGER;
    unsigned int ad = hash_combine( djb2(ID), hash(32*p + c) );
    while(index[ad]) {
      if( pos[index[ad] - 1] == p && chr[index[ad] -1] == c 
)//          && string_equal<T, SEXP>(ID, id[index[ad] - 1])) 
        return index[ad];
      ++ad %= m;
    }
    return NA_INTEGER;
  }

  // look up by snpid:chr:pos:alleles
  template<typename T1, typename T2>
  inline unsigned int lookup(T1 ID, int c, int p, T2 AL1, T2 AL2) const {
    if(htype != snpid_chr_pos_al)
      return NA_INTEGER;
    unsigned int ad = hash_combine( djb2(ID), hash(32*p + c) );
    while(index[ad]) {
      if( pos[index[ad] - 1] == p && chr[index[ad] -1] == c 
          && string_equal(ID, id[index[ad] - 1]) 
          && string_equal(AL1, A1[index[ad] - 1])
          && string_equal(AL2, A2[index[ad] - 1]))
        return index[ad];
      ++ad %= m;
    }
    return NA_INTEGER;
  }

  // look up by chr:pos
  inline unsigned int lookup(int c, int p) const {
    unsigned int ad = hash(32*p + c) ;
    while(index[ad]) {
      if(pos[index[ad] - 1] == p && chr[index[ad] -1] == c)
        return index[ad];
      ++ad %= m;
    }
    return NA_INTEGER;
  }

  // look up by chr:pos:alleles 
  template<typename T>
  inline unsigned int lookup(int c, int p, T AL1, T AL2) const {
    if(htype != chr_pos_al && htype != snpid_chr_pos_al)
      return NA_INTEGER;
    unsigned int ad = hash(32*p + c) ;
    while(index[ad]) {
      if( pos[index[ad] - 1] == p && chr[index[ad] -1] == c 
          && string_equal(AL1, A1[index[ad] - 1])
          && string_equal(AL2, A2[index[ad] - 1]))
        return index[ad];
      ++ad %= m;
    }
    return NA_INTEGER;
  }

  // look up by chr:pos:alleles 
  // trying strand flip or allele swap
  template<typename T>
  inline unsigned int lookup(int c, int p, T AL1, T AL2, bool & flip, bool & swap) const {
    if(htype != chr_pos_al && htype != snpid_chr_pos_al)
      return NA_INTEGER;
    unsigned int ad = hash(32*p + c) ;
    while(index[ad]) {
      if( pos[index[ad] - 1] == p && chr[index[ad] -1] == c) {
        if(string_equal(AL1, A1[index[ad] - 1]) && string_equal(AL2, A2[index[ad] - 1])) {
          flip = false; 
          swap = false;
          return index[ad];
        }
        std::string AL1_f = flip_strand(AL1);
        std::string AL2_f = flip_strand(AL2);
        if(string_equal(AL1_f, A1[index[ad] - 1]) && string_equal(AL2_f, A2[index[ad] - 1])) {
          flip = true; 
          swap = false;
          return index[ad];
        }
        if(string_equal(AL2, A1[index[ad] - 1]) && string_equal(AL1, A2[index[ad] - 1])) {
          flip = false; 
          swap = true;
          return index[ad];
        }
        if(string_equal(AL2_f, A1[index[ad] - 1]) && string_equal(AL1_f, A2[index[ad] - 1])) {
          flip = true; 
          swap = true;
          return index[ad];
        }
      }
      ++ad %= m;
    }
    return NA_INTEGER;
  }



  // ---------------- hash functions -------------------------
  inline unsigned int hash(unsigned int x) const {
    return (3141592653U * ((unsigned int)(x)) >> (32 - k));
  }

  inline unsigned int djb2(const char *str) const {
    unsigned int h = 5381;
    int c;
    while ((c = (unsigned char) *str++))
      h = (h << 5) + h + c; /*  33h + c */
    // Rcout << "djb2 value for string " << str << " = " << h << "\n";
    return hash(h);
  }

  inline unsigned int djb2(const std::string & str) const {
    return djb2(str.c_str());
  }

  inline unsigned int djb2(const SEXP str) const {
    return djb2(CHAR(str));
  }

  inline unsigned int hash_combine(unsigned int x, unsigned y) const {
    // Rcout << "combine x = " << x << " and y = " << y << " -> ";
    int h = x^y;
    // Rcout << "resultat " << h << "\n";
    return h;
  }

};
#endif

