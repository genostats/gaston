#include <RcppEigen.h>
#include <cmath>
#include <iostream>
#include "matrix-varia.h"
#ifndef GASTONLOGIT
#define GASTONLOGIT
using namespace Rcpp;
using namespace Eigen;

template<typename T1, typename T2>
void logistic_model(const Eigen::MatrixBase<T1> & y, const Eigen::MatrixBase<T2> & x, double eps, VectorXd & beta, MatrixXd & XWX_i) {
  int n(y.rows()), p(x.cols());
  VectorXd W(n);
  MatrixXd XWX(p,p), WX(n,p), P(n,n);
  VectorXd pi(n), z(n), Pz(n);
  double d, log_d;

  W.setZero();
  VectorXd U(p);
  double U_norm(1);
  beta.setZero();
  while (U_norm > eps) {
    for(int j = 0; j < n; j++) {
      pi(j) = 1/( 1 + exp( - x.row(j).dot(beta) ) );
      W(j) = pi(j)*(1-pi(j));
    }
    U.noalias() = x.transpose()*(y - pi);

    WX.noalias() = W.asDiagonal() * x;
    XWX.noalias() = x.transpose() * WX;
    sym_inverse(XWX, XWX_i, log_d, d, 1e-5);

    if(std::abs(d) < 1e-5) {
      for(int i = 0; i < p; i++) {
        beta(i) = NAN;
        for(int j = 0; j < p; j++) XWX_i(i,j) = NAN;
      }
      return;
    }

    beta += XWX_i*U;
    U_norm = U.norm();
  } 
}

// en float
template<typename T1, typename T2>
void logistic_model_f(const Eigen::MatrixBase<T1> & y, const Eigen::MatrixBase<T2> & x, float eps, VectorXf & beta, MatrixXf & XWX_i) {
  int n(y.rows()), p(x.cols());
  VectorXf W(n);
  MatrixXf XWX(p,p), WX(n,p), P(n,n);
  VectorXf pi(n), z(n), Pz(n);
  float d, log_d;
  W.setZero();
  VectorXf U(p);
  float U_norm(1);
  beta.setZero();
  int k = 0;
  while(U_norm > eps) {
    for(int j = 0; j < n; j++) {
      pi(j) = 1/( 1 + exp( - x.row(j).dot(beta) ) );
      W(j) = pi(j)*(1-pi(j));
    }
    U.noalias() = x.transpose()*(y - pi);

    WX.noalias() = W.asDiagonal() * x;
    XWX.noalias() = x.transpose() * WX;
    sym_inverse(XWX, XWX_i, log_d, d, 1e-5);

    if(std::abs(d) < 1e-5) {
      for(int i = 0; i < p; i++) {
        beta(i) = NAN;
        for(int j = 0; j < p; j++) XWX_i(i,j) = NAN;
      }
      return;
    }

    beta += XWX_i*U;
    U_norm = U.norm();
    if(k++ > 10) return;
  }
}
#endif
