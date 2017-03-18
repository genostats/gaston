#include <RcppEigen.h>
#include <iostream>
#include "matrix-varia.h"
using namespace Rcpp;
using namespace Eigen;

// p = nb de PC à mettre en effet fixe
// P0y n'est pas un vrai paramètre, c'est un output [les n-p dernières composantes de P0y]
// idem pour v, XViX_i [ correspond à (X_b' V_{0b}^{-1} X_b)^{-1} ]
template<typename T1, typename T2, typename T3>
void diago_likelihood_newton(double & h2, double & v, int p, const Eigen::MatrixBase<T1> & Y, Eigen::MatrixBase<T2> & X, const Eigen::MatrixBase<T3> & Sigma, 
                        VectorXd & P0y, MatrixXd & XViX_i, bool constraint, bool verbose) {
  int n = Sigma.rows();
  int r = X.cols();
  double log_d, d;
  double yP0y, y_PDP_y, y_PDPDP_y;
  
  VectorXd Deltab = Sigma.bottomRows(n-p) - VectorXd::Ones(n-p);

  for(int i = 0; i < 10; i++) {
    VectorXd V0b = h2*Sigma.bottomRows(n-p) + (1-h2)*VectorXd::Ones(n-p); 
    VectorXd V0bi = V0b.cwiseInverse();
 
    MatrixXd ViX = V0bi.asDiagonal() * X.bottomRows(n-p);  //      V_{0b}^{-1} X_b
    MatrixXd XViX = X.bottomRows(n-p).transpose() * ViX;   // X_b' V_{0b}^{-1} X_b
    XViX_i = MatrixXd(r,r);
    sym_inverse(XViX, XViX_i, log_d, d, 1e-5);

    // P0y = P0b * Y.bottomRows(n-p);
    P0y = V0bi.asDiagonal() * Y.bottomRows(n-p) - ViX * ( XViX_i * ( ViX.transpose() * Y.bottomRows(n-p) ) ) ;

    yP0y = P0y.dot( Y.bottomRows(n-p).col(0) );
    
    VectorXd DeltaP0y = Deltab.cwiseProduct( P0y );
    y_PDP_y = P0y.dot( DeltaP0y );
  
    // VectorXd P0DeltaP0y = P0b * DeltaP0y;
    VectorXd P0DeltaP0y  = V0bi.asDiagonal() * DeltaP0y - ViX * ( XViX_i * ( ViX.transpose() * DeltaP0y ) ) ;
    y_PDPDP_y = DeltaP0y.dot( P0DeltaP0y );


    // double f = -0.5*(V0b.array().log().sum() + log_d + (n-r-p)*log(yP0y)+ (n-r-p)*(1-log((double)(n-r-p))));
    // **** calcul dérivées
    // MatrixXd P0b = V0bi.asDiagonal();
    // P0b.noalias() -= ViX * XViX_i * ViX.transpose();
    // MatrixXd DeltaP0 = Deltab.asDiagonal() * P0b;
    // double df0 = DeltaP0.trace() - (n-r-p)*y_PDP_y/yP0y;
    // double ddf0 = -trace_of_product(DeltaP0, DeltaP0) + (n-r-p)*( 2*y_PDPDP_y/yP0y - (y_PDP_y*y_PDP_y)/(yP0y*yP0y));

    // alternative qui évite le calcul de P0b et Delta P0b (qui sont rapides à calculer
    // mais de grande taille ->O(n²r) remplacé par O(r²n) si je ne m'abuse
    MatrixXd B = ViX.transpose() * Deltab.asDiagonal() * ViX;
    MatrixXd C = ViX.transpose() * (Deltab.cwiseProduct(V0bi.cwiseProduct(Deltab))).asDiagonal() * ViX;
    MatrixXd AB = XViX_i * B;
    double trace_PD = Deltab.cwiseProduct(V0bi).sum() - AB.trace();
    double trace_PDPD = Deltab.cwiseProduct(V0bi.cwiseProduct(Deltab.cwiseProduct(V0bi))).sum() - 2*trace_of_product(XViX_i, C) + trace_of_product(AB,AB);

    double df = trace_PD - (n-r-p)*y_PDP_y/yP0y;
    double ddf = -trace_PDPD + (n-r-p)*( 2*y_PDPDP_y/yP0y - (y_PDP_y*y_PDP_y)/(yP0y*yP0y));
    if(verbose) {
      Rcout << "[Iteration " << i+1 << "] ";
      Rcout << "h2 = " << h2 << " df = " << -0.5*df;
    }

    h2 -= df/ddf;
  }
  v = yP0y / (n-r-p);

}

