#include <iostream>
#include <Rcpp.h>
#include <stdint.h>
#define DEBUG false

using namespace Rcpp;

// matrices qui stockent de façon compacte des valeurs de 1 à 4
class proxy_matrix4 {
  public:
    proxy_matrix4(uint8_t * r, size_t b) : ref(r), j(b) {}
    // set value
    inline int operator=(int x) { 
      // Rcout << "proxy set, j = " << j << ", int x = " << x << "\n";
      uint8_t & a = ref[j/4];
      a &= ~(3 << ((j%4)*2));  // set to 00
      a |= (x << ((j%4)*2)); // set to x
      return x;
    }
    inline int operator=(proxy_matrix4 x) {
      // Rcout << "proxy set, j = " << j << "(x is a proxy)\n";
      *this = (int) x;
      return (int) x;
    }
  
    // get value
    inline operator int() const {
      // Rcout << "proxy (int), j = " << j << "\n";
      return((int) ((ref[j/4] >> ((j%4)*2)) & 3));
    }
  private:  
    uint8_t * ref;
    size_t j;
};


class matrix4 {
  public:
    // constructor
    matrix4(size_t a = 0, size_t b = 0);
    // copy constructor
    matrix4(const matrix4&);
    // copy constructor from numeric matrix
    matrix4(const NumericMatrix);
    matrix4(const RawMatrix);
    // destructor
    ~matrix4();
    // affectation
    matrix4& operator=(const matrix4&);
    matrix4& operator=(const NumericMatrix);

    // set and get
    inline uint8_t get(size_t i, size_t j) const;
    // void set(size_t i, size_t j, int val);
    inline void set(size_t i, size_t j, uint8_t val);

    // ()
    inline uint8_t operator()(size_t i, size_t j) const;
    proxy_matrix4 operator()(size_t i, size_t j) {
      // Rcout << "cree proxy i = " << i << ", j = " << j << "\n";
      return proxy_matrix4(data[i], j);
    }

    void fill_line(size_t li, NumericVector w);

    // output 
    friend std::ostream& operator<<(std::ostream&, const matrix4);

    size_t nrow, ncol;
    size_t true_ncol;
    uint8_t ** data;
    void allocations() {
      data = new uint8_t * [nrow];
      for(size_t i = 0; i < nrow; i++) {
        data[i] = new uint8_t [true_ncol];
        // on initialise avec des 3 partout (3 -> NA)
        // important pour accélérer certaines fonctions que la matrice soit bordée avec des NA
        std::fill(data[i], data[i]+true_ncol, 255);
      }
    }
};

// set and get
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

