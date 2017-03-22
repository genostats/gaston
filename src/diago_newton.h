#include <RcppEigen.h>
#include <iostream>
#include "matrix-varia.h"
#include "optimize_2.h"

using namespace Rcpp;
using namespace Eigen;

// structure qui hérite fun pour l'optimisation par Brent
// Normalement l'utilisation doit être
// T_Matrix = MatrixXd, T_Vector = VectorXd, T_scal = double
// T_Matrix = MatrixXf, T_Vector = VectorXf, T_scal = float
template<typename T_Matrix, typename T_Vector, typename T_scal>
class diag_likelihood : public fun {
  public:
    int p, n, r;
    T_Matrix Y;
    T_Matrix X;
    T_Matrix Sigma;
    T_Vector P0y;
    T_scal & v;
    T_Matrix XViX_i;
    T_Vector Deltab;
    T_scal likelihood;
    T_scal d, log_d;
    T_Vector V0b, V0bi;
    T_Matrix ViX, XViX;
    T_scal yP0y;
    diag_likelihood(int p, const T_Matrix & Y, const T_Matrix & X, const T_Vector & Sigma, 
                   T_Vector & P0y, T_scal & v, T_Matrix & XViX_i) : p(p), n(Sigma.rows()), r(X.cols()), Y(Y), X(X), Sigma(Sigma), P0y(P0y), v(v), XViX_i(XViX_i) {
      Deltab = Sigma.bottomRows(n-p) - T_Vector::Ones(n-p);
    } ;

    void update(T_scal h2) {
      V0b = h2*Sigma.bottomRows(n-p) + (1-h2)*T_Vector::Ones(n-p); 
      V0bi = V0b.cwiseInverse();
      ViX = V0bi.asDiagonal() * X.bottomRows(n-p);  //      V_{0b}^{-1} X_b
      XViX = X.bottomRows(n-p).transpose() * ViX;   // X_b' V_{0b}^{-1} X_b
      sym_inverse(XViX, XViX_i, log_d, d, 1e-5);

      P0y = V0bi.asDiagonal() * Y.bottomRows(n-p) - ViX * ( XViX_i * ( ViX.transpose() * Y.bottomRows(n-p) ) ) ;
      yP0y = P0y.dot( Y.bottomRows(n-p).col(0) );
      v = yP0y / (n-r-p);
    }

    T_scal f(T_scal h2) {
      update(h2);
      likelihood = -0.5*(V0b.array().log().sum() + log_d + (n-r-p)*log(yP0y) + (n-r-p)*(1-log((T_scal)(n-r-p))) ); // avec le terme constant
      return likelihood;
    }

    void df_ddf(T_scal h2, T_scal & df, T_scal & ddf) {
      update(h2);
      T_Vector DeltaP0y = Deltab.cwiseProduct( P0y );
      T_scal y_PDP_y = P0y.dot( DeltaP0y );
  
      // T_Vector P0DeltaP0y = P0b * DeltaP0y;
      T_Vector P0DeltaP0y  = V0bi.asDiagonal() * DeltaP0y - ViX * ( XViX_i * ( ViX.transpose() * DeltaP0y ) ) ;
      T_scal y_PDPDP_y = DeltaP0y.dot( P0DeltaP0y );

      // T_scal f = -0.5*(V0b.array().log().sum() + log_d + (n-r-p)*log(yP0y)+ (n-r-p)*(1-log((T_scal)(n-r-p))));
      // **** calcul dérivées
      // T_Matrix P0b = V0bi.asDiagonal();
      // P0b.noalias() -= ViX * XViX_i * ViX.transpose();
      // T_Matrix DeltaP0 = Deltab.asDiagonal() * P0b;
      // T_scal df0 = DeltaP0.trace() - (n-r-p)*y_PDP_y/yP0y;
      // T_scal ddf0 = -trace_of_product(DeltaP0, DeltaP0) + (n-r-p)*( 2*y_PDPDP_y/yP0y - (y_PDP_y*y_PDP_y)/(yP0y*yP0y));

      // alternative qui évite le calcul de P0b et Delta P0b (qui sont rapides à calculer
      // mais de grande taille ->O(n²r) remplacé par O(r²n) si je ne m'abuse
      T_Matrix B = ViX.transpose() * Deltab.asDiagonal() * ViX;
      T_Matrix C = ViX.transpose() * (Deltab.cwiseProduct(V0bi.cwiseProduct(Deltab))).asDiagonal() * ViX;
      T_Matrix AB = XViX_i * B;
      T_scal trace_PD = Deltab.cwiseProduct(V0bi).sum() - AB.trace();
      T_scal trace_PDPD = Deltab.cwiseProduct(V0bi.cwiseProduct(Deltab.cwiseProduct(V0bi))).sum() - 2*trace_of_product(XViX_i, C) + trace_of_product(AB,AB);
  
      df = trace_PD - (n-r-p)*y_PDP_y/yP0y;
      ddf = -trace_PDPD + (n-r-p)*( 2*y_PDPDP_y/yP0y - (y_PDP_y*y_PDP_y)/(yP0y*yP0y));
    }
    void newton(T_scal & h2, const T_scal min_h2, const T_scal max_h2, const T_scal eps, const bool verbose);
};

template<typename T_Matrix, typename T_Vector, typename T_scal>
void diag_likelihood<T_Matrix, T_Vector, T_scal>::newton(T_scal & h2, const T_scal min_h2, const T_scal max_h2, const T_scal eps, const bool verbose) {
  int nb_reseeds = 0, i = 0;
  T_scal df = 1+2*eps; 

  bool tried_max = false, tried_min = false;
  if(h2 == min_h2) tried_min = true;
  if(h2 == max_h2) tried_max = true;

  while(std::abs(df) > 2*eps) {

    T_scal ddf;
    df_ddf(h2, df, ddf);

    if(verbose) {
      Rcout << "[Iteration " << ++i << "] ";
      Rcout << "h² = " << h2 << " df = " << -0.5*df << std::endl;
    }

    // si on est au bord de l'intervalle et qu'on ne tend pas à revenir dedans
    if(h2 == min_h2 && !isnan(df) && df > 0) {
      if(verbose) Rcout << "[Iteration " << i << "] maximum at min h² = " << h2 << std::endl;
      break;
    }
    if(h2 == max_h2 && !isnan(df) && df < 0) {
      if(verbose) Rcout << "[Iteration " << i << "] maximum at max h² = " << h2 << std::endl;
      break;
    }

    // si la convexité est mauvaise
    if(ddf < 0) {
      if(verbose) Rcout << "[Iteration " << i << "] likelihood isn't concave" << std::endl;
      T_scal old_h2 = h2;
      if(df < 0) {
        if(!tried_max) {
          h2 = max_h2; 
          tried_max = true;
          if(verbose) Rcout << "[Iteration " << i << "] restarting from h2 = " << h2 << std::endl;
          continue;
        } 
        if(verbose) Rcout << "[Iteration " << i << "] Using Brent algorithm" << std::endl;
        h2 = Brent_fmin(old_h2, max_h2, 1e-5);
        // pour mettre à jour P0Y etc... [est-ce bien utile ? le dernier point où est calculé f doit suffire]
        update(h2); 
        if(verbose) Rcout << "[Iteration " << i << "] Brent gives h² = " << h2 << std::endl;
        break;
      }
      if(df > 0) { 
        if(!tried_min) {
          h2 = min_h2;
          tried_min = true;
          if(verbose) Rcout << "[Iteration " << i << "] restarting from h2 = " << h2 << std::endl;
          continue;
        }
        if(verbose) Rcout << "[Iteration " << i << "] Using Brent algorithm" << std::endl;
        Brent_fmin(min_h2, old_h2, 1e-5);
        // pour mettre à jour P0Y etc...
        update(h2);
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
        return;
      }
    } else if(h2 < min_h2) {
      h2 = min_h2;
      tried_min = true;
    } else if(h2 > max_h2) {
      h2 = max_h2;
      tried_max = true;
    }
  }
  Rcout << "A.P0y size " << P0y.size() << "\n";
  v = yP0y / (n-r-p);
  Rcout << "ok then\n";
}

