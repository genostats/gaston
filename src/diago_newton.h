#include <RcppEigen.h>
#include <iostream>
#include "matrix-varia.h"
#include "optimize_2.h"

using namespace Rcpp;
using namespace Eigen;

// structure qui hérite fun pour l'optimisation par Brent
// Normalement l'utilisation doit être
// MATRIX = MatrixXd, VECTOR = VectorXd, scalar_t = double
// MATRIX = MatrixXf, VECTOR = VectorXf, scalar_t = float
template<typename MATRIX, typename VECTOR, typename scalar_t>
class diag_likelihood : public fun<scalar_t> {
  public:
    int p, n, r;
    MATRIX Y;
    MATRIX X;
    MATRIX Sigma;
    VECTOR P0y;
    scalar_t & v;
    MATRIX XViX_i;
    VECTOR Deltab;
    scalar_t likelihood;
    scalar_t d, log_d;
    VECTOR V0b, V0bi;
    MATRIX ViX, XViX;
    scalar_t yP0y;
    diag_likelihood(int p, const MATRIX & Y, const MATRIX & X, const VECTOR & Sigma, 
                   VECTOR & P0y, scalar_t & v, MATRIX & XViX_i) : p(p), n(Sigma.rows()), r(X.cols()), Y(Y), X(X), Sigma(Sigma), P0y(P0y), v(v), XViX_i(XViX_i) {
      Deltab = Sigma.bottomRows(n-p) - VECTOR::Ones(n-p);
    } ;

    void update(scalar_t h2) {
      V0b = h2*Sigma.bottomRows(n-p) + (1-h2)*VECTOR::Ones(n-p); 
      V0bi = V0b.cwiseInverse();
      ViX = V0bi.asDiagonal() * X.bottomRows(n-p);  //      V_{0b}^{-1} X_b
      XViX = X.bottomRows(n-p).transpose() * ViX;   // X_b' V_{0b}^{-1} X_b
      sym_inverse(XViX, XViX_i, log_d, d, 1e-5);

      P0y = V0bi.asDiagonal() * Y.bottomRows(n-p) - ViX * ( XViX_i * ( ViX.transpose() * Y.bottomRows(n-p) ) ) ;
      yP0y = P0y.dot( Y.bottomRows(n-p).col(0) );
      v = yP0y / (n-r-p);
    }

    scalar_t f(scalar_t h2) {
      update(h2);
      likelihood = -0.5*(V0b.array().log().sum() + log_d + (n-r-p)*log(yP0y) + (n-r-p)*(1-log((scalar_t)(n-r-p))) ); // avec le terme constant
      return likelihood;
    }

    void df_ddf(scalar_t h2, scalar_t & df, scalar_t & ddf) {
      update(h2);
      VECTOR DeltaP0y = Deltab.cwiseProduct( P0y );
      scalar_t y_PDP_y = P0y.dot( DeltaP0y );
  
      // VECTOR P0DeltaP0y = P0b * DeltaP0y;
      VECTOR P0DeltaP0y  = V0bi.asDiagonal() * DeltaP0y - ViX * ( XViX_i * ( ViX.transpose() * DeltaP0y ) ) ;
      scalar_t y_PDPDP_y = DeltaP0y.dot( P0DeltaP0y );

      // scalar_t f = -0.5*(V0b.array().log().sum() + log_d + (n-r-p)*log(yP0y)+ (n-r-p)*(1-log((scalar_t)(n-r-p))));
      // **** calcul dérivées
      // MATRIX P0b = V0bi.asDiagonal();
      // P0b.noalias() -= ViX * XViX_i * ViX.transpose();
      // MATRIX DeltaP0 = Deltab.asDiagonal() * P0b;
      // scalar_t df0 = DeltaP0.trace() - (n-r-p)*y_PDP_y/yP0y;
      // scalar_t ddf0 = -trace_of_product(DeltaP0, DeltaP0) + (n-r-p)*( 2*y_PDPDP_y/yP0y - (y_PDP_y*y_PDP_y)/(yP0y*yP0y));

      // alternative qui évite le calcul de P0b et Delta P0b (qui sont rapides à calculer
      // mais de grande taille ->O(n²r) remplacé par O(r²n) si je ne m'abuse
      MATRIX B = ViX.transpose() * Deltab.asDiagonal() * ViX;
      MATRIX C = ViX.transpose() * (Deltab.cwiseProduct(V0bi.cwiseProduct(Deltab))).asDiagonal() * ViX;
      MATRIX AB = XViX_i * B;
      scalar_t trace_PD = Deltab.cwiseProduct(V0bi).sum() - AB.trace();
      scalar_t trace_PDPD = Deltab.cwiseProduct(V0bi.cwiseProduct(Deltab.cwiseProduct(V0bi))).sum() - 2*trace_of_product(XViX_i, C) + trace_of_product(AB,AB);
  
      df = trace_PD - (n-r-p)*y_PDP_y/yP0y;
      ddf = -trace_PDPD + (n-r-p)*( 2*y_PDPDP_y/yP0y - (y_PDP_y*y_PDP_y)/(yP0y*yP0y));
    }

    void beta_blup(scalar_t h2, VECTOR & beta, VECTOR & omega) {
      update(h2);
      beta_blup(beta, omega);
    }

    // *********** CALCUL DES BLUPS ************************
    // Attention P0y n'est que (P0y)b, les n-p dernières composantes ! (les p premières sont nulles)
    void beta_blup(VECTOR & beta, VECTOR & omega) {
      VECTOR Sigmab = sigma.bottomRows(n-p);
      VECTOR omega = h2 * sigmab.asDiagonal() * P0y;

      VECTOR z = y;
      z.tail(n-p) -= omega + (1-h2)*P0y;
      // Xb' Xb
      MATRIX xtx( MATRIX(r,r).setZero().selfadjointView<Lower>().rankUpdate( x.bottomRows(n-p).transpose() ));
      MATRIX xtx0( xtx );
      MATRIX xtxi(r,r); // et son inverse
      double d1, ld1;
      sym_inverse(xtx0, xtxi, d1, ld1, 1e-5); // détruit xtx0

      VECTOR beta(r+p);
      beta.topRows(r) = xtxi * x.bottomRows(n-p).transpose() * z.bottomRows(n-p);
      beta.bottomRows(p) = z.topRows(p) - x.topRows(p) * beta.topRows(r);
    }

    void newton(scalar_t & h2, const scalar_t min_h2, const scalar_t max_h2, const scalar_t eps, const bool verbose);
};




template<typename MATRIX, typename VECTOR, typename scalar_t>
void diag_likelihood<MATRIX, VECTOR, scalar_t>::newton(scalar_t & h2, const scalar_t min_h2, const scalar_t max_h2, const scalar_t eps, const bool verbose) {
  int nb_reseeds = 0, i = 0;
  scalar_t df = 1+2*eps; 

  bool tried_max = false, tried_min = false;
  if(h2 == min_h2) tried_min = true;
  if(h2 == max_h2) tried_max = true;

  while(std::abs(df) > 2*eps) {

    scalar_t ddf;
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
      scalar_t old_h2 = h2;
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

