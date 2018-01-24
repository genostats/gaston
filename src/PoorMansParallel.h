#include <thread>
#include <vector>
#include <Rcpp.h>

#ifndef POOR_MANS_PARALLEL
#define POOR_MANS_PARALLEL


namespace PoorMansParallel {

  /********* From RcppParallel/Common.h *************/
  struct Worker {  
     // construct and destruct (delete virtually)
     Worker() {}
     virtual ~Worker() {}
   
     // dispatch work over a range of values
     virtual void operator()(std::size_t begin, std::size_t end) = 0;   
   
     // disable copying and assignment
  private:
     Worker(const Worker&);
     void operator=(const Worker&);
  };
  // Tag type used for disambiguating splitting constructors
  struct Split {};

  /********* end From **/

  // private elements
  namespace {
    // number of threads, with default initial value
    unsigned int n = std::thread::hardware_concurrency();

    // this is useful because std::thread can't handle a
    // link to a member function
    template <typename Reducer>
    void run_operator(Reducer * r, size_t a, size_t b) {
      r->operator()(a,b);
    }
  } // fin private

  // access to n (defined in PoorMansParallel.cpp)
  unsigned int nb_threads();
  void set_nb_threads(unsigned int nb);
  void set_nb_threads_to_default();

  // split function for class that inherit Worker
  template <typename Reducer>
  inline Reducer * split(Reducer & A) {
    Split s;
    return new Reducer(A, s);
  }

  // *************************************************************** //
  // Poor Man's Parallel Reduce
  template <typename Reducer>
  void parallelReduce(size_t begin, size_t end, Reducer& A) {
    if(end <= begin) // nothing to do
      return;  
   unsigned int nb_threads = PoorMansParallel::nb_threads();
   size_t le = end-begin;
   size_t size = le/nb_threads;
  
   std::vector<Reducer*> red(nb_threads-1);
   std::vector<std::thread> thr(nb_threads-1);
 
   // split A and start thread 
   size_t be = begin;
   for(size_t i = 0; i < nb_threads-1; i++) {
     size_t en = be + size;
     red[i] = split(A);
     thr[i] = std::thread( run_operator<Reducer>, red[i], be, en);                                       
     be = en;
   }
   // don't forget last block
   A(be, end);
 
   // join threads and Reducers, and destruct splitted Reducers
   for(size_t i = 0; i < nb_threads-1; i++) {
     thr[i].join();
     A.join(*red[i]);
     delete red[i];
   }
 }
  // Poor Man's Parallel For (defined in .cpp)
  void parallelFor(size_t begin, size_t end, Worker& A);
  void parallelFor(size_t begin, size_t end, Worker& A, int chunk);

  // chunk is ignored
  template <typename Reducer>
  void parallelReduce(size_t begin, size_t end, Reducer& A, int chunk) {
    parallelReduce<Reducer>(begin, end, A);
  }

} // end namespace
 
#endif
