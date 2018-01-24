#include "PoorMansParallel.h"

namespace PoorMansParallel {
  // access to n
  unsigned int nb_threads() {
    return n;
  }

  void set_nb_threads(unsigned int nb) {
    n = nb;
  }

  void set_nb_threads_to_default() {
    n = std::thread::hardware_concurrency();
  }

  // Poor Man's Parallel For
  void parallelFor(size_t begin, size_t end, Worker& A) {
    if(end <= begin) // nothing to do
      return;
   unsigned int nb_threads = PoorMansParallel::nb_threads();
   size_t le = end-begin;
   size_t size = le/nb_threads;

   std::vector<std::thread> thr(nb_threads-1);

   // start threads
   size_t be = begin;
   for(size_t i = 0; i < nb_threads-1; i++) {
     size_t en = be + size;
     thr[i] = std::thread( run_operator<Worker>, &A, be, en);
     be = en;
   }
   // don't forget last block
   A(be, end);

   // join threads
   for(size_t i = 0; i < nb_threads-1; i++) {
     thr[i].join();
   }
 }
  // chunk is ignored
  void parallelFor(size_t begin, size_t end, Worker& A, int chunk) {
    parallelFor(begin, end, A);
  }
} // end namespace



// functions to export to R...

int _get_nb_threads() {
  // return std::thread::hardware_concurrency();
  return PoorMansParallel::nb_threads();
}

void _set_nb_threads(int n) {
  PoorMansParallel::set_nb_threads(n);
}

//[[Rcpp::export]]
void _set_nb_threads_to_default() {
  PoorMansParallel::set_nb_threads_to_default();
}
 
RcppExport SEXP get_nb_threads() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(_get_nb_threads());
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP set_nb_threads(SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    _set_nb_threads(n);
    return R_NilValue;
END_RCPP
}

RcppExport SEXP set_nb_threads_to_default() {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    _set_nb_threads_to_default();
    return R_NilValue;
END_RCPP
}


