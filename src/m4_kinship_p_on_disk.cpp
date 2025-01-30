#include <Rcpp.h>
#include "Parallel.h"
#include <iostream>
#include "matrix4.h"
#include "m4_kinship_type.h"
#include "mmatrix.h"
#include "mmatrix_methods.h"


#include <cstdlib> // for rand() in creating K workers
#include <cstdio> //for deleting files
#include <cstring> //for opening errno if error

using namespace Rcpp;
using namespace Parallel;


// ************* [on ne symmétrise pas] ***********

struct paraKin_p : public Worker {
  // input and others
  uint8_t ** data;
  const size_t ncol;
  const size_t true_ncol;
  const Ktype n;                  // should be (nb_snps - 1)
  const double * p;
  const size_t sizeK;
  const bool dominance;
  MMatrix<Ktype> *K_file;

  // output
  Ktype * K;
  
  // constructeurs
  paraKin_p(uint8_t ** data, const size_t ncol, const size_t true_ncol, const Ktype n, const double * p, const bool dominance) :
          data(data), ncol(ncol), true_ncol(true_ncol), n(n), p(p), sizeK((4*true_ncol)*(4*true_ncol+1)/2) , dominance(dominance) { 
          //generate a random name for K
          std::cout << "sizeK = " << sizeK << std::endl;
          std::string path = "paraKin_worker1_" + std::to_string(rand());
          // TO DEBUG
          std::cout << "creating a K with 1st constructor, called " << path << "\n";
          K_file = new MMatrix<Ktype>(path, sizeK, 1); // K is padded to a multiple of 4...
          K = K_file->data();
          //no need for std::fill, mmatrix already full of zeroes
        }
  paraKin_p(paraKin_p & Q, Split) : data(Q.data), ncol(Q.ncol), true_ncol(Q.true_ncol), n(Q.n), p(Q.p), sizeK(Q.sizeK), dominance(Q.dominance) {
          //TODO JU : add handling of K on disk (un pointeur vers une array de Ktype dynamically allocated)
          //generate a random name for K
          std::cout << "sizeK = " << sizeK << std::endl;
          std::string path = "paraKin_worker2_" + std::to_string(rand());
          // TO DEBUG
          std::cout << "creating a K with 2nd constructor, called " << path << "\n";
          K_file = new MMatrix<Ktype>(path, sizeK, 1);
          K = K_file->data();
        }

  // destructeur
  // TODO JU : remove the K files also !
  ~paraKin_p() { 
    std::cout << "Trying to DELETE a paraKIn_p\n";
    if (K_file) {
      std::cout << "K_file exists, checking path...\n";
      try {
        std::string  path = K_file->path();
        std::cout << "Deleting file: " << path << "\n";
      } catch (...) {
            std::cerr << "Error accessing K_file->path()\n";
      }
      //std::cout << "Deleting file: " << path << "\n";
      delete K_file;// will flush before unmapping
      //int result = std::remove(path.c_str());// so deleting file after matrix
      //if (result != 0) printf( "%s\n", strerror( errno ) );
    }
  }

  // worker !
  void operator()(size_t beg, size_t end) {
    Ktype H[4];
    Ktype h0[32];
    Ktype h1[32];
    H[3] = 0;

    for(size_t i = beg; i < end; i++) {
      Ktype p_ = (Ktype) p[i];   // freq de l'allele 1 (noté q le plus souvent !!)
      if(p_ == 0 || p_ == 1) continue;
      Ktype v0, v1, v2;
      if(dominance) {
        v0 = p_/(1-p_);
        v1 = -1;
        v2 = (1-p_)/p_; 
      } else {
        Ktype w_ = 1/sqrt(2*p_*(1-p_)); 
        v0 = -2*p_*w_;
        v1 = (1-2*p_)*w_;
        v2 = (2-2*p_)*w_;
      }
      // on divise par le nombre de SNPs
      H[0] = v0/sqrt(n);
      H[1] = v1/sqrt(n);
      H[2] = v2/sqrt(n);
      
      size_t k = 0;
      for(size_t a1 = 0; a1 < 4; a1++) {
        for(size_t a2 = 0; a2 < 4; a2++) {
          h0[k++] = H[a2];
          h0[k++] = H[a1];
        }
      }

      uint8_t * dd = data[i];
      
      k = 0;
      for(size_t j1 = 0; j1 < true_ncol; j1++) {
        uint8_t x1 = dd[j1];
        for(unsigned int ss1 = 0; (ss1 < 4); ss1++) {
          for(size_t a = 0; a < 32; a++) h1[a] = H[x1&3]*h0[a];
          for(size_t j2 = 0; j2 < j1; j2++) {
            uint8_t x2 = dd[j2];
            K[k++] += h1[ (x2&15)<<1 ];
            K[k++] += h1[ ((x2&15)<<1)+1 ];
            K[k++] += h1[ (x2>>4)<<1 ];
            K[k++] += h1[ ((x2>>4)<<1)+1 ];
          }
          size_t j2 = j1;
          uint8_t x2 = dd[j2];
          for(unsigned int ss2 = 0; ss2 <= ss1; ss2++) {
            K[k++] += H[x1&3]*H[x2&3]; // h1[(x2&3)<<1];
            x2 >>= 2;
          } 
          x1 >>= 2; 
        }
      } 
    }
  }

  // recoller
  void join(const paraKin_p & Q) {
    std::transform(K, K + sizeK, Q.K, K, std::plus<Ktype>());
    // autrement dit : K += Q.K;
  }

};

//[[Rcpp::export]]
std::string Kinship_pw_on_disk(XPtr<matrix4> p_A, const std::vector<double> & p, LogicalVector snps, bool dominance, int chunk) {
  // Will return the name of the file containing the numeric matrix with results

    int nb_snps = sum(snps);

    if (snps.length() != p_A->nrow || p.size() != nb_snps)
        stop("Dimensions mismatch");

    uint8_t **data = new uint8_t *[nb_snps];
    size_t k = 0;
    for (size_t i = 0; i < p_A->nrow; i++)
    {
        if (snps[i])
            data[k++] = p_A->data[i];
    }

  // TODO JU : add handling of K in x on disk (generate new_names)
  paraKin_p X(data, p_A->ncol, p_A->true_ncol, (Ktype) (nb_snps - 1), &p[0], dominance);
  parallelReduce(0, nb_snps, X, chunk);

  delete [] data;

  // hardcoded double type because of the casting in the for loop after
  MMatrix<double> Y_disk("Kinship_matrix", p_A->ncol, p_A->ncol);

  k = 0;
  for(size_t i = 0; i < p_A->ncol; i++) {
    for(size_t j = 0; j <= i; j++) {
      Y_disk(j,i) = (double) X.K[k++];
    }
  }

  // symmetriser
  k = 0;
  for(size_t i = 0; i < p_A->ncol; i++) {
    for(size_t j = 0; j <= i; j++) {
      Y_disk(i,j) = (double) X.K[k++]; // ou Y(j,i)
    }
  }

  return "Kinship_matrix";
}

void populate_bed_matrix_JU(XPtr<matrix4> ref) {
    uint8_t new_data = 1;

    for (size_t i = 0; i < ref->nrow; i++) {
        for (size_t j = 0; j < ref->ncol; j++) {
            ref->set(i, j, new_data);
            new_data = (new_data + 1) % 256;
            // TO DEBUVG
            //uint8_t k = ref->data[i][j];
            //printf("value got from ref->data at indices i = %zu, j = %zu , %u\n",i,j,k);
        }
    }
}

void head_bed_matrix_JU(XPtr<matrix4> ref) {
  if (ref->nrow > 9 && ref->ncol > 9)
  {
    for (size_t j = 0; j < 10; j++)
    {
      printf("\nrow %zu : ", j);
      for (size_t i = 0; i < 10; i++)
      {
        printf("%u,  ", ref->get(i,j));
      }
    }
  }
  else {
    printf("Warning : printing head on less than 10x10 matrix ! \n");
    for (size_t j = 0; j < ref->nrow; j++)
    {
      printf("\nrow %zu : ", j);
      for (size_t i = 0; i < ref->ncol; i++)
      {
        printf("%u,  ", ref->get(i,j));
      }
    }
  }
}


void copy_bed_mat_on_disk_JU(XPtr<matrix4> ref) {
  
  // TODO JU : to modify with this value :
  //true_ncol = ncol/4 + ((ncol%4 == 0u)?0:1);
  
  //populate_bed_matrix_JU(ref);

  MMatrix<uint8_t> Test_disk_uint8_t("Test_disk", ref->ncol, ref->nrow);

  printf("\n\n----------uint8_t's loop:---------------\n\n");
  uint8_t k = 0;

  for (size_t j = 0; j < ref->nrow; j++)
  {
    for (size_t i = 0; i < ref->ncol; i++)
    {
      k = ref->get(j, i);
      Test_disk_uint8_t.at(j, i) = k;
    }
  }
  printf("\n\n----------After the loop (in ref):---------------\n\n");
  head_bed_matrix_JU(ref);
  }