#include <iostream>
#include <Rcpp.h>
#ifndef GASTONloubar
#define GASTONloubar

using namespace Rcpp;

/* Classes de vecteurs et de matrices très rustiques
 * pour éviter RcppArmadillo qui ralentit le code parallèle sur le cluster
 * on ne met pas grd chose dedans (intérêt #1 : avoir un destructeur
 * pour ne pas avoir à gérer ça dans les fonctions)
 */

#ifndef DEBUG 
#define DEBUG false
#endif

struct D {
  D() {}
} ; // désambiguation constructeur
const D Clone;


// une classe de vecteur très rustique
class bar {
  public:
  bar() : n(0), mine(false)  { 
    #if DEBUG
    Rcout << "bar()\n";  
    #endif
  }

  bar(size_t n) : n(n), mine(true) { 
    #if DEBUG
    Rcout << "bar(" << n << ")\n";  
    #endif
    data = new double [n];
 }

  bar(const bar & A) : n(A.n), mine(true) { 
    #if DEBUG
    Rcout << "bar(bar & A) (n = " << n << ") (copie)\n";
    #endif
    data = new double[n];
    for(size_t i = 0; i < n; i++) data[i] = A.data[i];
  }
  
  bar(NumericVector A) : n(A.length()), mine(true) {
    #if DEBUG
    Rcout << "bar(NumericVector) length = " << n << "(copie)\n";
    #endif
    data = new double[n];
    double *pA = &A[0];
    for(size_t i = 0; i < n; i++) data[i] = pA[i];
  }

  bar & operator=(const bar & A) { 
    #if DEBUG
    Rcout << "bar = (copie)\n";
    #endif
    if(A.n != n) {
      this->~bar();
      n = A.n;
      data = new double [n];
      mine = true;
    }
    for(size_t i = 0; i < n; i++) data[i] = A.data[i];
    return *this;
  }

  bar(NumericVector A, D) : n(A.length()), mine(false), data(&A[0]) {
    #if DEBUG
    Rcout << "bar(NumericVector length = " << n << ", Clone)\n";
    #endif
  }

  bar(double * A, size_t n) : n(n), mine(true) {
    #if DEBUG
    Rcout << "bar(double *, n = " << n << ") (copie)\n";
    #endif
    data = new double [n];
    for(size_t i = 0; i < n; i++) data[i] = A[i];
  }

  bar(double * A, size_t n, D) : n(n), mine(false), data(A) {
    #if DEBUG
    Rcout << "bar(double *, n = " << n << ", Clone)\n";
    #endif
  }

  ~bar() { 
    if(mine) { 
      #if DEBUG
      Rcout << "~bar mine(" << n << ")\n"; 
      #endif
      delete [] data; 
    } else { 
      #if DEBUG
      Rcout << "~bar not mine(" << n << ")\n";
      #endif
    } 
  }

  inline void zeros() { std::fill(data, data + n, 0); }
 
  size_t n;
  bool mine;
  double * data;
};

// classe de matrices construite sur la classe de vecteurs
// c'est super rustique, c'est à l'utilisateur d'accéder
// à l'élément (i,j) en faisant (par exemple) data[i+j*ncol]
class lou : public bar {
  public: 
  lou(size_t nrow, size_t ncol) : nrow(nrow), ncol(ncol), n(nrow*ncol), x(bar(nrow*ncol)), data(x.data) {}
  lou(const lou & a) : nrow(a.nrow), ncol(a.ncol), n(nrow*ncol), x(a.x), data(x.data) {}
  lou(NumericMatrix A) : nrow(A.nrow()), ncol(A.ncol()), n(nrow*ncol), x(A), data(x.data) {}
  lou(NumericMatrix A, D d) : nrow(A.nrow()), ncol(A.ncol()), n(nrow*ncol), x(A, d), data(x.data) {}
  lou(double * A, size_t nrow, size_t ncol) : nrow(nrow), ncol(ncol),  n(nrow*ncol),x(A, nrow*ncol), data(x.data) {}
  lou(double * A, size_t nrow, size_t ncol, D d) : nrow(nrow), ncol(ncol), n(nrow*ncol), x(A, nrow*ncol, d), data(x.data) {}
  
  // le destructeur par défaut fait l'affaire...
  size_t nrow;
  size_t ncol;
  size_t n;
  bar x;
  double * data;
};

// classe de ragged array
// tjs très rustique.
// élément (i,j) -> .data[csl[i]+j]
class rag : public bar {
  public:
  rag() : nrow(0), n(0), mine(false) {}

  rag(size_t nrow, size_t * le_) : nrow(nrow), mine(true) {
    #if DEBUG 
    Rcout << "rag(" << nrow << ",...)\n";
    #endif
    if(nrow > 0) {
      le = new size_t [nrow];
      cs = new size_t [nrow];
      cs[0] = 0;
      le[0] = le_[0];
      #if DEBUG
      Rcout << "rag() le[0] = " << le[0] << "\n";
      #endif
      for(size_t i = 1; i < nrow; i++) {
        le[i] = le_[i];  
        #if DEBUG
        Rcout << "rag() le[i] = " << le[i] << "\n";
        #endif
        cs[i] = cs[i-1] + le[i-1];
      }
      n = cs[nrow-1]+le[nrow-1];
      data = new double [n];
    }
  }

  // un constructeur où l'utilisateur est prié de 
  // mettre à jour le reste plus tard
  rag(size_t nrow) : nrow(nrow), mine(false) {
    #if DEBUG
    Rcout << "rag(" << nrow << ")\n";
    #endif
    if(nrow > 0) {
      le = new size_t [nrow];
      cs = new size_t [nrow];
    }
  }

  void set_length(size_t * le_) {
   if(nrow > 0) {
      cs[0] = 0;
      le[0] = le_[0];
      #if DEBUG
      Rcout << "rag() le[0] = " << le[0] << "\n";
      #endif
      for(size_t i = 1; i < nrow; i++) {
        le[i] = le_[i];  
        // Rcout << "rag() le[i] = " << le[i] << "\n";
        cs[i] = cs[i-1] + le[i-1];
      }
      mine = true;
      n = cs[nrow-1]+le[nrow-1];
      data = new double [n];
    }
  }

  ~rag() {
    #if DEBUG
    Rcout << "~rag()\n";
    #endif
    if(nrow > 0) {
      delete [] cs;
      delete [] le;
    }
    if(mine) delete [] data;
  }

  size_t nrow;
  size_t * le;  // length
  size_t * cs; // cum sum length

  size_t n;
  double * data;
  bool mine;
};

/******** produit lou x bar ****/

void loubar(const lou & A, const bar & X, bar & R);
bar loubar(const lou & A, const bar & X);
void barlou(const bar & X, const lou & A, bar & R);
bar barlou(const bar & X, const lou & A);

// bar x bar (produit scalaire)
double barbar(const bar & X, const bar & Y);

// lou x lou
void loulou(const lou & A, const lou & B, lou & C);

/* conversion vers objets R */

NumericMatrix as_r(lou & A);

NumericVector as_r(bar & A);


#endif
