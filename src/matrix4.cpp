#include <Rcpp.h>
#include <iostream>
#include "matrix4.h"

using namespace Rcpp;

// constructeur
matrix4::matrix4(size_t a, size_t b) : nrow(a), ncol(b) {
  #if DEBUG
  Rcout << "constructeur (" << a << ", " << b << ")\n";
  #endif
  true_ncol = b/4 + ((b%4 == 0)?0:1);
  allocations();
}

// constructeur par copie
matrix4::matrix4(const matrix4& x) : nrow(x.nrow), ncol(x.ncol), true_ncol(x.true_ncol) {
  #if DEBUG
  Rcout << "constructeur par copie (" << x.nrow << ", " << x.ncol << ")\n";
  #endif
  allocations();
  for(size_t i = 0; i < nrow; i++) {
    for(size_t j = 0; j < true_ncol; j++) data[i][j] = x.data[i][j];
  }
}


// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// LES COPIES DES MATRICES R FONT UNE TRANSPOSITION
matrix4::matrix4(NumericMatrix x) {
  #if DEBUG
  Rcout << "constructeur par copie Numeric Matrix (" << x.nrow() << ", " << x.ncol() << ")\n";
  #endif
  ncol = x.nrow();
  nrow = x.ncol();
  true_ncol = ncol/4 + ((ncol%4 == 0u)?0:1);
  allocations();
  for(size_t i = 0; i < nrow; i++) {
    for(size_t j = 0; j < ncol; j++) {
      int u = NumericMatrix::is_na(x(j,i))?3:(int) x(j,i);
      set(i,j, (u>>2 == 0?u:3));
    }
  }
} 

matrix4::matrix4(RawMatrix x) {
  #if DEBUG
  Rcout << "constructeur par copie RawMatrix (" << x.nrow() << ", " << x.ncol() << ")\n";
  #endif
  ncol = x.nrow();
  nrow = x.ncol();
  true_ncol = ncol/4 + ((ncol%4 == 0u)?0:1);
  allocations();
  for(size_t i = 0; i < nrow; i++) {
    for(size_t j = 0; j < ncol; j++) {
      uint8_t u = NumericMatrix::is_na(x(j,i))?3:(uint8_t) x(j,i);
      set(i,j, (u>>2 == 0?u:3));
    }
  }
} 
  
// destructeur
matrix4::~matrix4() { 
  #if DEBUG
  Rcout << "destruction nrow = " << nrow << ", ncol = " << ncol << "\n";
  #endif
  for(size_t i = 0; i < nrow; i++) delete[] data[i];
  if(nrow > 0) delete[] data;
}

// affectation
matrix4& matrix4::operator=(const matrix4& x) {
  #if DEBUG
  Rcout << "copie (affectation)\n";
  #endif
  if(&x != this) {
    // si les dimensions diffèrent, détruire et recréer...
    if(nrow != x.nrow || true_ncol != x.true_ncol) {
      this->~matrix4(); 
      nrow = x.nrow;
      true_ncol = x.true_ncol;
      allocations();
    }
    ncol = x.ncol; 
    for(size_t i = 0; i < nrow; i++) {
      for(size_t j = 0; j < true_ncol; j++) data[i][j] = x.data[i][j];
    }
  }
  return *this;
}

// copie d'une NumericMatrix
matrix4& matrix4::operator=(const NumericMatrix x) {
  if(nrow != (size_t) x.nrow() || ncol != (size_t) x.ncol()) {
    this->~matrix4();
    nrow = x.nrow();
    ncol = x.ncol();
    true_ncol = ncol/4 + ((ncol%4 == 0u)?0:1);
    allocations();
  }
  for(size_t i = 0; i < nrow; i++) {
    for(size_t j = 0; j < ncol; j++) {
      int u = NumericMatrix::is_na(x(i,j))?3:(int) x(i,j);
      this->set(i,j, (u>>2 == 0?u:3));
    }
  }
  return *this;
}

// pour la lecture d'un fichier ligne à ligne
// (attention, sans doute nécessaire de transposer...)
// Ici contrairement à set on teste les valeurs (cette fonction sera donnée à l'utilisateur final)
void matrix4::fill_line(size_t li, NumericVector w) {
  if((size_t) w.length() != ncol) {
    Rcout << "fill_line : Length mismatch, nothing done\n";
    return;
  }
  if(li >= nrow) {
    Rcout << "fill_line : Line number " << li << "too high (should be between 0 and " << nrow-1 << ")\n";
    return;
  }
  // on efface la ligne
  std::fill(data[li], data[li]+true_ncol, 255);
  for(size_t i = 0; i < true_ncol - 1; i++) {
    uint8_t &a = data[li][i];      
    for(int ss = 0; ss < 4; ss++) {
      a <<= 2;
      uint8_t x = NumericVector::is_na(w[4*i+3-ss])?3:(uint8_t) w[4*i+3-ss];
      a |= (x >> 2 == 0?x:3);
    }
  }
  size_t i = true_ncol - 1;
  uint8_t &a = data[li][i];
  for(int ss = 4*i+4-ncol; ss < 4; ss++) {
    a <<= 2;
    uint8_t x = NumericVector::is_na(w[4*i+3-ss])?3:(uint8_t) w[4*i+3-ss];
    // uint8_t x = (uint8_t) w(4*i+3-ss);
    a |= ((x >> 2 == 0)?x:3);
  }
}


// set and get 
/*
uint8_t matrix4::get(size_t i, size_t j) const {
  return((int) ((data[i][j/4] >> ((j%4)*2)) & 3));
}

void matrix4::set(size_t i, size_t j, uint8_t val) {
  uint8_t & a = data[i][j/4];
  a &= ~(3 << ((j%4)*2));  // set to 00
  a |= (val << ((j%4)*2)); // set to val
}

uint8_t matrix4::operator()(size_t i, size_t j) const {
  #if DEBUG
  Rcout << "(const) int(), i = " << i << ", j = " << j << "n";
  #endif
  return get(i,j);
}
*/

// <<
std::ostream& operator<<(std::ostream& o, const matrix4 x) {
  for(size_t i = 0; i < x.nrow; i++) {
    o << "[" << i << ",] ";
    for(size_t j = 0; j < x.ncol; j++) {
      o << (int) x.get(i,j);
    }
    o << "\n";
  }
  return o;
}

