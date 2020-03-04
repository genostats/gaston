#include <RcppEigen.h>
#include <cmath>
#include <iostream>
#ifndef GASTONLOGITMODEL
#define GASTONLOGITMODEL

template<typename scalar_t>
using MATRIX = Eigen::Matrix<scalar_t, Eigen::Dynamic, Eigen::Dynamic>;

template<typename scalar_t>
using VECTOR = Eigen::Matrix<scalar_t, Eigen::Dynamic, 1>;

using namespace Rcpp;
using namespace Eigen;

// on ne peut pas utiliser a.log().sum() parce qu'avec scalar_t = float ça donne un warning de compilation.
// -> workaround....
template<typename scalar_t>
inline scalar_t sum_of_log(const VECTOR<scalar_t> & a) {
  scalar_t S(0);
  int n = a.rows();
  for(Eigen::Index i = 0; i < n; ++i) S += std::log(a(i));
  return S;
}

template<typename scalar_t>
void logistic_model2(const VECTOR<scalar_t> & y, const MATRIX<scalar_t> & x, VECTOR<scalar_t> & beta, MATRIX<scalar_t> & XWX_i, 
                     scalar_t eps = 1.0e-8, int max_iter = 25) {
  int n(y.rows()), p(x.cols());
  VECTOR<scalar_t> W(n), pi(n), U(p), Xbeta(n);
  MATRIX<scalar_t> XWX(p,p), WX(n,p);
  scalar_t dev(0), dev_old(0);

  W.setZero();

  beta.setZero();
  for(int k = 0; k < max_iter; k++) {
    Xbeta.noalias() = x * beta;
    for(int i = 0; i < n; i++) pi(i) = 1/( 1 + exp( - Xbeta(i) ) );
    W.noalias() = pi.cwiseProduct(VECTOR<scalar_t>::Ones(n) - pi);

    U.noalias() = x.transpose()*(y - pi);

    WX.noalias() = W.asDiagonal() * x;
    XWX.noalias() = x.transpose() * WX;

    beta += XWX.llt().solve(U);  // better than sym_inverse(XWX, XWX_i, log_d, d, 1e-5);

    // deviance
    dev_old = dev;
    // dev = y.dot(Xbeta) + (VECTOR<scalar_t>::Ones(n) - pi).log().sum();
    dev = y.dot(Xbeta) + sum_of_log<scalar_t>(VECTOR<scalar_t>::Ones(n) - pi);

    if(k > 0 && fabs(dev - dev_old)/(fabs(dev) + 0.1) < eps) 
      break; 
  } 
  // Mais il nous faut l'inverse de XWX pour calculer l'écart type
  XWX_i = XWX.llt().solve( MATRIX<scalar_t>::Identity(p,p) );
}


// THE EXACT SAME CODE EXCEPT FOR THE OFFSET ...
template<typename scalar_t>
void logistic_model_offset(const VECTOR<scalar_t> & y, const VECTOR<scalar_t> & offset, const MATRIX<scalar_t> & x, VECTOR<scalar_t> & beta, 
                           MATRIX<scalar_t> & XWX_i, scalar_t eps = 1.0e-8, int max_iter = 25) {
  int n(y.rows()), p(x.cols());
  VECTOR<scalar_t> W(n), pi(n), U(p), Xbeta(n);
  MATRIX<scalar_t> XWX(p,p), WX(n,p);
  scalar_t dev(0), dev_old(0);

  W.setZero();

  beta.setZero();
  for(int k = 0; k < max_iter; k++) {
    Xbeta.noalias() = x * beta + offset;
    for(int i = 0; i < n; i++) pi(i) = 1/( 1 + exp( - Xbeta(i) ) );
    W.noalias() = pi.cwiseProduct(VECTOR<scalar_t>::Ones(n) - pi);

    U.noalias() = x.transpose()*(y - pi);

    WX.noalias() = W.asDiagonal() * x;
    XWX.noalias() = x.transpose() * WX;

    beta += XWX.llt().solve(U);  // better than sym_inverse(XWX, XWX_i, log_d, d, 1e-5);

    // deviance
    dev_old = dev;
    dev = y.dot(Xbeta) + (VECTOR<scalar_t>::Ones(n) - pi).log().sum();

    if(k > 0 && fabs(dev - dev_old)/(fabs(dev) + 0.1) < eps) 
      break; 
  } 
  // Mais il nous faut l'inverse de XWX pour calculer l'écart type
  XWX_i = XWX.llt().solve( MATRIX<scalar_t>::Identity(p,p) );
}
#endif
