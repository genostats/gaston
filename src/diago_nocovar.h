#include <RcppEigen.h>
#include <iostream>
#include "matrix-varia.h"
using namespace Rcpp;
using namespace Eigen;

// p = nb de PC à mettre en effet fixe
// P0y n'est pas un vrai paramètre, c'est un output [les n-p dernières composantes de P0y]
// ICI, cas sans covariable
template<typename T1, typename T3>
double diago_likelihood(double h2, int p, const Eigen::MatrixBase<T1> & Y, const Eigen::MatrixBase<T3> & Sigma, 
                        VectorXd & P0y, double & v) {
  int n = Sigma.rows();
  double yP0y;

  VectorXd V0b = h2*Sigma.bottomRows(n-p) + (1-h2)*VectorXd::Ones(n-p); 
  VectorXd V0bi = V0b.cwiseInverse();

  // Py = P0y / v [seulement les n-p dernieres compostantes sont calculées, les p premières sont nulles]
  // Ici on a simple V0b^{-1} dans le coin sud est...
  P0y = V0bi.asDiagonal() * Y.bottomRows(n-p);
  yP0y = P0y.dot( Y.bottomRows(n-p).col(0) );
  v = yP0y / (n-p);
  Rcout << "h2 = " << h2 << ", v = " << v << "\n";

  return -0.5*(V0b.array().log().sum() + (n-p)*log(yP0y) + (n-p)*(1-log((double)(n-p))) );
}

template<typename T1, typename T3>
double diago_likelihood(double tau, double s2, int p, Eigen::MatrixBase<T1> & Y,
                        Eigen::MatrixBase<T3> & Sigma, VectorXd & P0y) {
  int n = Sigma.rows();
  double yP0y;
  
  double v  = s2 + tau;
  double h2 = tau/v;

  VectorXd V0b = h2*Sigma.bottomRows(n-p) + (1-h2)*VectorXd::Ones(n-p); 
  VectorXd V0bi = V0b.cwiseInverse();

  P0y = V0bi.asDiagonal() * Y.bottomRows(n-p);
  yP0y = P0y.dot( Y.bottomRows(n-p).col(0) );

  return -0.5*(V0b.array().log().sum() + yP0y/v + (n-p)*log(v));
}



