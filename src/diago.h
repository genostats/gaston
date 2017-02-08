#include <RcppEigen.h>
#include <iostream>
#include "matrix-varia.h"
using namespace Rcpp;
using namespace Eigen;

// p = nb de PC à mettre en effet fixe
// P0y n'est pas un vrai paramètre, c'est un output [les n-p dernières composantes de P0y]
// idem pour v, XViX_i [ correspond à (X_b' V_{0b}^{-1} X_b)^{-1} ]
template<typename T1, typename T2, typename T3>
double diago_likelihood(double h2, int p, const Eigen::MatrixBase<T1> & Y, Eigen::MatrixBase<T2> & X, const Eigen::MatrixBase<T3> & Sigma, 
                        VectorXd & P0y, double & v, MatrixXd & XViX_i) {
  int n = Sigma.rows();
  int r = X.cols();
  double log_d, d;
  double yP0y;

  VectorXd V0b = h2*Sigma.bottomRows(n-p) + (1-h2)*VectorXd::Ones(n-p); 
  VectorXd V0bi = V0b.cwiseInverse();

  MatrixXd ViX = V0bi.asDiagonal() * X.bottomRows(n-p);  //      V_{0b}^{-1} X_b
  MatrixXd XViX = X.bottomRows(n-p).transpose() * ViX;   // X_b' V_{0b}^{-1} X_b
  XViX_i = MatrixXd(r,r);
  sym_inverse(XViX, XViX_i, log_d, d, 1e-5);

  // Py = P0y / v [seulement les n-p dernieres compostantes sont calculées, les p premières sont nulles]
  P0y = V0bi.asDiagonal() * Y.bottomRows(n-p) - ViX * XViX_i * ViX.transpose() * Y.bottomRows(n-p);
  yP0y = P0y.dot( Y.bottomRows(n-p).col(0) );
  v = yP0y / (n-r-p);

  return -0.5*(V0b.array().log().sum() + log_d + (n-r-p)*log(yP0y) + (n-r-p)*(1-log((double)(n-r-p))) ); // on met le terme constant pour pouvoir comparer à celle ci-dessous
}

template<typename T1, typename T2, typename T3>
double diago_likelihood(double tau, double s2, int p, Eigen::MatrixBase<T1> & Y, Eigen::MatrixBase<T2> & X, 
                        Eigen::MatrixBase<T3> & Sigma, VectorXd & P0y, MatrixXd & XViX_i) {
  int n = Sigma.rows();
  int r = X.cols();
  double log_d, d;
  double yP0y;
  
  double v  = s2 + tau;
  double h2 = tau/v;

  VectorXd V0b = h2*Sigma.bottomRows(n-p) + (1-h2)*VectorXd::Ones(n-p); 
  VectorXd V0bi = V0b.cwiseInverse();

  MatrixXd ViX = V0bi.asDiagonal() * X.bottomRows(n-p);  //      V_{0b}^{-1} X_b
  MatrixXd XViX = X.bottomRows(n-p).transpose() * ViX;   // X_b' V_{0b}^{-1} X_b
  XViX_i = MatrixXd(r,r);
  sym_inverse(XViX, XViX_i, log_d, d, 1e-5);

  // Py = P0y / v [seulement les n-p dernieres compostantes sont calculées, les p premières sont nulles]
  P0y = V0bi.asDiagonal() * Y.bottomRows(n-p) - ViX * XViX_i * ViX.transpose() * Y.bottomRows(n-p);
  yP0y = P0y.dot( Y.bottomRows(n-p).col(0) );

  return -0.5*(V0b.array().log().sum() + log_d + yP0y/v + (n-r-p)*log(v));
}




