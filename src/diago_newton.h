#include <RcppEigen.h>
#include <iostream>
#include "matrix-varia.h"
#include "optimize_2.h"

using namespace Rcpp;
using namespace Eigen;

// structure contenant tous les paramètres pour diago_likelihood
template<typename T1, typename T2, typename T3>
class diag_likelihood : public fun {
    int p;
    MatrixXd Y;
    MatrixXd X;
    MatrixXd Sigma;
    VectorXd P0y;
    double v;
    MatrixXd XViX_i;
    double likelihood;
  public:
    diag_likelihood(int p, const Eigen::MatrixBase<T1> & Y, const Eigen::MatrixBase<T2> & X, const Eigen::MatrixBase<T3> & Sigma, 
                   VectorXd & P0y, double & v, MatrixXd & XViX_i) : p(p), Y(Y), X(X), Sigma(Sigma), P0y(P0y), v(v), XViX_i(XViX_i) {} ;
    double f(double h2) {
      int n = Sigma.rows();
      int r = X.cols();
      double log_d, d;
      double yP0y;
      Rcout << "h2 = " << h2;

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
      likelihood = -0.5*(V0b.array().log().sum() + log_d + (n-r-p)*log(yP0y) + (n-r-p)*(1-log((double)(n-r-p))) ); // avec le terme constant
      Rcout << " likelihood = " << likelihood << "\n";
      return likelihood;
    }
};


template<typename T1, typename T2, typename T3>
void inline diago_likelihood_derivatives(double & df, double & ddf, double & h2, int n, int r, int p, VectorXd & Deltab, 
                                         const Eigen::MatrixBase<T1> & Y,
                                         const Eigen::MatrixBase<T2> & X, const Eigen::MatrixBase<T3> & Sigma, VectorXd & P0y, 
                                         MatrixXd & XViX_i, double & yP0y) {
    double log_d, d;
    VectorXd V0b = h2*Sigma.bottomRows(n-p) + (1-h2)*VectorXd::Ones(n-p); 
    VectorXd V0bi = V0b.cwiseInverse();
 
    MatrixXd ViX = V0bi.asDiagonal() * X.bottomRows(n-p);  //      V_{0b}^{-1} X_b
    MatrixXd XViX = X.bottomRows(n-p).transpose() * ViX;   // X_b' V_{0b}^{-1} X_b
    sym_inverse(XViX, XViX_i, log_d, d, 1e-5);

    // P0y = P0b * Y.bottomRows(n-p);
    P0y = V0bi.asDiagonal() * Y.bottomRows(n-p) - ViX * ( XViX_i * ( ViX.transpose() * Y.bottomRows(n-p) ) ) ;

    yP0y = P0y.dot( Y.bottomRows(n-p).col(0) );
    
    VectorXd DeltaP0y = Deltab.cwiseProduct( P0y );
    double y_PDP_y = P0y.dot( DeltaP0y );
  
    // VectorXd P0DeltaP0y = P0b * DeltaP0y;
    VectorXd P0DeltaP0y  = V0bi.asDiagonal() * DeltaP0y - ViX * ( XViX_i * ( ViX.transpose() * DeltaP0y ) ) ;
    double y_PDPDP_y = DeltaP0y.dot( P0DeltaP0y );

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

    df = trace_PD - (n-r-p)*y_PDP_y/yP0y;
    ddf = -trace_PDPD + (n-r-p)*( 2*y_PDPDP_y/yP0y - (y_PDP_y*y_PDP_y)/(yP0y*yP0y));
}

// p = nb de PC à mettre en effet fixe
// P0y n'est pas un vrai paramètre, c'est un output [les n-p dernières composantes de P0y]
// idem pour v, XViX_i [ correspond à (X_b' V_{0b}^{-1} X_b)^{-1} ]
template<typename T1, typename T2, typename T3>
int diago_likelihood_newton(double & h2, double & v, int p, const Eigen::MatrixBase<T1> & Y, const Eigen::MatrixBase<T2> & X, const Eigen::MatrixBase<T3> & Sigma, 
                            VectorXd & P0y, MatrixXd & XViX_i, const double min_h2, const double max_h2, const double eps, const bool verbose) {
  int n = Sigma.rows();
  int r = X.cols();
  int nb_reseeds = 0, i = 0;
  double yP0y;
  double df = 1+2*eps; 
  Rcout << "min = " << min_h2 << " max = " << max_h2 << "\n";

  VectorXd Deltab = Sigma.bottomRows(n-p) - VectorXd::Ones(n-p);
  XViX_i = MatrixXd(r,r);

  bool tried_max = false, tried_min = false;
  while(std::abs(df) > 2*eps) {

    if(verbose) {
      Rcout << "[Iteration " << ++i << "] ";
      Rcout << "h² = " << h2 << " df = " << -0.5*df << std::endl;
    }
    double ddf;

    diago_likelihood_derivatives(df, ddf, h2, n, r, p, Deltab, Y, X, Sigma, P0y, XViX_i, yP0y);

    // si on est au bord de l'intervalle et qu'on ne tend pas à revenir dedans
    if(h2 == min_h2 && !isnan(df) && df > 0) {
      if(verbose) Rcout << "[Iteration " << i << "] maximum at min h² = " << h2 << std::endl;
      break;
    }
    if(h2 == max_h2 && !isnan(df) && df < 0) {
      if(verbose) Rcout << "[Iteration " << i << "] maximum at max h² = " << h2 << " " << max_h2 << std::endl;
      break;
    }

    // si la convexité est mauvaise
    if(ddf < 0) {
      if(verbose) Rcout << "[Iteration " << i << "] likelihood isn't concave" << std::endl;
      double old_h2 = h2;
      if(df < 0) {
        if(!tried_max) {
          h2 = max_h2;
          tried_max = true;
          diago_likelihood_derivatives(df, ddf, h2, n, r, p, Deltab, Y, X, Sigma, P0y, XViX_i, yP0y);
          if(df < 0) {
            if(verbose) Rcout << "[Iteration " << i << "] maximum at max h² = " << h2 << std::endl;
            break; 
          } else if(ddf > 0) {
            if(verbose) Rcout << "[Iteration " << i << "] restarting from h2 = " << h2 << std::endl;
            continue;
          }
        } 
        if(verbose) Rcout << "[Iteration " << i << "] Using Brent algorithm" << std::endl;
        diag_likelihood< T1, T2, T3 > A(p, Y, X, Sigma, P0y, v, XViX_i);
        h2 = A.Brent_fmin(old_h2, max_h2, 1e-5);
        diago_likelihood_derivatives(df, ddf, h2, n, r, p, Deltab, Y, X, Sigma, P0y, XViX_i, yP0y);
        if(verbose) Rcout << "[Iteration " << i << "] Brent gives h² = " << h2 << std::endl;
        break;
      }
      if(df > 0) { 
        if(!tried_min) {
          h2 = min_h2;
          tried_min = true;
          diago_likelihood_derivatives(df, ddf, h2, n, r, p, Deltab, Y, X, Sigma, P0y, XViX_i, yP0y);
          if(df > 0) {
            if(verbose) Rcout << "[Iteration " << i << "] maximum at h² = " << h2 << std::endl;
            break; 
          } else if(ddf > 0) {
            if(verbose) Rcout << "[Iteration " << i << "] restarting from h2 = " << h2 << std::endl;
            continue;
          }
        }
        if(verbose) Rcout << "[Iteration " << i << "] Using Brent algorithm" << std::endl;
        diag_likelihood< T1, T2, T3 > A(p, Y, X, Sigma, P0y, v, XViX_i);
        h2 = A.Brent_fmin(min_h2, old_h2, 1e-5);
        diago_likelihood_derivatives(df, ddf, h2, n, r, p, Deltab, Y, X, Sigma, P0y, XViX_i, yP0y);
        if(verbose) Rcout << "[Iteration " << i << "] Brent gives h² = " << h2 << std::endl;
        break;
      }
    }
    h2 -= df/ddf;

    if(std::isnan(h2)) {
      if(nb_reseeds++ < 5) {
        h2 = R::runif(min_h2,max_h2);
        if(verbose) Rcout << "[Iteration " << i << "] restarting from random value h² = " << h2 << std::endl;
      } else {
        if(verbose) Rcout << "[Iteration " << i << "] canceling optimization" << std::endl;
        return 0;
      }
    } else if(h2 < min_h2) {
      h2 = min_h2;
      tried_min = true;
    } else if(h2 > max_h2) {
      h2 = max_h2;
      tried_max = true;
    }
  }

  /*if(constraint) {
    bool recompute = false;
    if(h2 < 0) { 
      h2 = 0; recompute = true;
      if(verbose) Rcout << "Constraining h² = 0" << std::endl;
    }
    if(h2 > 1) {
      h2 = 1; recompute = true;
      if(verbose) Rcout << "Constraining h² = 1" << std::endl;
    }
    if(recompute) { // il faut recalculer XViX_i P0y et YP0y 
      VectorXd V0b = h2*Sigma.bottomRows(n-p) + (1-h2)*VectorXd::Ones(n-p); 
      VectorXd V0bi = V0b.cwiseInverse();
      MatrixXd ViX = V0bi.asDiagonal() * X.bottomRows(n-p);  //      V_{0b}^{-1} X_b
      MatrixXd XViX = X.bottomRows(n-p).transpose() * ViX;   // X_b' V_{0b}^{-1} X_b
      double log_d, d;
      sym_inverse(XViX, XViX_i, log_d, d, 1e-5);

      P0y = V0bi.asDiagonal() * Y.bottomRows(n-p) - ViX * ( XViX_i * ( ViX.transpose() * Y.bottomRows(n-p) ) ) ;
      yP0y = P0y.dot( Y.bottomRows(n-p).col(0) );
    }
  } */
  v = yP0y / (n-r-p);
  return 1;
}

