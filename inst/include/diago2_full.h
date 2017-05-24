#include <RcppEigen.h>
#include <iostream>
#include "matrix-varia.h"
#include "optimize_2.h"
#ifndef GASTONdiag_full_likelihood
#define GASTONdiag_full_likelihood

using namespace Rcpp;
using namespace Eigen;

// structure qui hérite fun pour l'optimisation par Brent
// Normalement l'utilisation doit être
// MATRIX = MatrixXd, VECTOR = VectorXd, scalar_t = double
// MATRIX = MatrixXf, VECTOR = VectorXf, scalar_t = float
template<typename MATRIX, typename VECTOR, typename scalar_t>
class diag_full_likelihood : public fun<scalar_t> {
  public:
    int p, n, r;
    const MATRIX Y;
    MATRIX X;
    const MATRIX Sigma;
    VECTOR P0y;
    scalar_t v;
    MATRIX XViX_i;
    VECTOR Delta, Deltab;
    scalar_t d, log_d;
    VECTOR V0, V0b, V0i, V0bi;
    MATRIX ViX, XViX, xtx;
    scalar_t yP0y;
    diag_full_likelihood(int p, const MATRIX & Y, const MATRIX & X, const VECTOR & Sigma) : 
       p(p), n(Sigma.rows()), r(X.cols()), Y(Y), X(X), Sigma(Sigma) {
         Delta = Sigma - VECTOR::Ones(n);
         Deltab = Delta.tail(n-p);
         XViX_i = MATRIX(r,r);
    } ;

    void update(scalar_t h2) {
      V0 =  h2*Sigma + (1-h2)*VECTOR::Ones(n);
      V0i = V0.cwiseInverse();
      V0b = V0.tail(n-p);
      V0bi = V0i.tail(n-p);
      ViX = V0bi.asDiagonal() * X.bottomRows(n-p);  //      V_{0b}^{-1} X_b
      XViX = X.bottomRows(n-p).transpose() * ViX;   // X_b' V_{0b}^{-1} X_b
      sym_inverse(XViX, XViX_i, log_d, d, 1e-5);

      P0y = V0bi.asDiagonal() * Y.bottomRows(n-p) - ViX * ( XViX_i * ( ViX.transpose() * Y.bottomRows(n-p) ) ) ;
      yP0y = P0y.dot( Y.bottomRows(n-p).col(0) );
      v = yP0y / n;
    }

    scalar_t f(scalar_t h2) {
      update(h2);
      return -0.5*(V0.array().log().sum() + n*log(yP0y) + n*(1-log((scalar_t) n)) ); // avec le terme constant
    }

    // vraisemblance au dernier point calculé
    scalar_t likelihood() {
      return -0.5*(V0.array().log().sum() + n*log(yP0y) + n*(1-log((scalar_t) n)) ); // avec le terme constant
    }

    scalar_t likelihood(scalar_t tau, scalar_t s2) {
      scalar_t vv = s2 + tau;
      scalar_t h2 = tau/vv;
      update(h2);
      return -0.5*(V0.array().log().sum() + yP0y/vv + (n)*log(vv));
    }
 
    void df_ddf(scalar_t h2, scalar_t & df, scalar_t & ddf) {
      update(h2);
      VECTOR DeltaP0y = Deltab.cwiseProduct( P0y );
      scalar_t y_PDP_y = P0y.dot( DeltaP0y );
  
      // VECTOR P0DeltaP0y = P0b * DeltaP0y;
      VECTOR P0DeltaP0y  = V0bi.asDiagonal() * DeltaP0y - ViX * ( XViX_i * ( ViX.transpose() * DeltaP0y ) ) ;
      scalar_t y_PDPDP_y = DeltaP0y.dot( P0DeltaP0y );

      VECTOR ViD = V0i.cwiseProduct(Delta);
      scalar_t trace_ViD = ViD.sum();
      scalar_t trace_ViDViD = ViD.cwiseProduct(ViD).sum();
  
      df = -0.5*(trace_ViD - n*y_PDP_y/yP0y);
      ddf = -0.5*(-trace_ViDViD + n*( 2*y_PDPDP_y/yP0y - (y_PDP_y*y_PDP_y)/(yP0y*yP0y)));
    }



    // *********** CALCUL DES BLUPS ************************
    // Attention P0y n'est que (P0y)b, les n-p dernières composantes ! (les p premières sont nulles)
    void blup(scalar_t h2, VECTOR & beta, VECTOR & omega, bool updateh2, bool only_first_r) {
      if(updateh2) update(h2);

      VECTOR Sigmab = Sigma.bottomRows(n-p);
      omega = h2 * Sigmab.asDiagonal() * P0y;

      VECTOR z = Y;
      z.tail(n-p) -= omega + (1-h2)*P0y;
      // Xb' Xb
      // MATRIX xtx( MATRIX(r,r).setZero().selfadjointView<Lower>().rankUpdate( X.bottomRows(n-p).transpose() ));
      xtx = MATRIX(r,r).setZero();
      SelfAdjointView<MATRIX, Lower> xtx_sa(xtx);
      xtx_sa.rankUpdate( X.bottomRows(n-p).transpose() );  
      MATRIX xtx0( xtx_sa );
      MATRIX xtxi(r,r); // et son inverse
      scalar_t d1, ld1;
      sym_inverse(xtx0, xtxi, d1, ld1, 1e-5); // détruit xtx0

      if(only_first_r) {
        beta = VECTOR(r);
        beta = xtxi * X.bottomRows(n-p).transpose() * z.bottomRows(n-p);
      } else {
        beta = VECTOR(r+p);
        beta.topRows(r) = xtxi * X.bottomRows(n-p).transpose() * z.bottomRows(n-p);
        beta.bottomRows(p) = z.topRows(p) - X.topRows(p) * beta.topRows(r);
     }
   }

};
#endif
