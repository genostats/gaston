#include <RcppEigen.h>

#ifndef GASTONANYNAN
#define GASTONANYNAN

template<typename scalar_t>
using VECTOR = Eigen::Matrix<scalar_t, Eigen::Dynamic, 1>;

template <typename T>
inline bool any_nan(VECTOR<T> x) {
  for(unsigned int i = 0; i < x.size(); i++) {
    if(std::isnan(x[i]))
      return true;
  }
  return false;
}
#endif
