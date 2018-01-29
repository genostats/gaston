#include <RcppEigen.h>
#include <iostream>
#include "matrix-varia.h"
#include "optimize_2.h"
#ifndef GASTONdiag_likelihood_nocovar
#define GASTONdiag_likelihood_nocovar
using namespace Rcpp;
using namespace Eigen;

// structure qui hérite fun pour l'optimisation par Brent
// Normalement l'utilisation doit être
// MATRIX = MatrixXd, VECTOR = VectorXd, scalar_t = double
// MATRIX = MatrixXf, VECTOR = VectorXf, scalar_t = float
template<typename MATRIX, typename VECTOR, typename scalar_t>
class diag_likelihood_nocovar : public fun<scalar_t> {
  public:
    int p, n;
    const MATRIX Y;
    const MATRIX Sigma;
    VECTOR P0y;
    scalar_t v;
    VECTOR Deltab;
    VECTOR V0b, V0bi;
    scalar_t yP0y;
    diag_likelihood_nocovar(int p, const MATRIX & Y, const VECTOR & Sigma) : 
       p(p), n(Sigma.rows()), Y(Y), Sigma(Sigma) {
         Deltab = Sigma.bottomRows(n-p) - VECTOR::Ones(n-p);
    } ;

    void update(scalar_t h2) {
      V0b = h2*Sigma.bottomRows(n-p) + (1-h2)*VECTOR::Ones(n-p); 
      V0bi = V0b.cwiseInverse();
      P0y = V0bi.asDiagonal() * Y.bottomRows(n-p);
      yP0y = P0y.dot( Y.bottomRows(n-p).col(0) );
      v = yP0y / (n-p);
    }

    scalar_t f(scalar_t h2) {
      update(h2);
      return -0.5*(V0b.array().log().sum() + (n-p)*log(yP0y) + (n-p)*(1-log((scalar_t)(n-p))) ); // avec le terme constant
    }
    
    // vraisemblance au dernier point calculé
    scalar_t likelihood() {
      return -0.5*(V0b.array().log().sum() + (n-p)*log(yP0y) + (n-p)*(1-log((scalar_t)(n-p))) ); // avec le terme constant
    }

    scalar_t likelihood(scalar_t tau, scalar_t s2) {
      scalar_t vv = s2 + tau;
      scalar_t h2 = tau/vv;
      update(h2);
      return -0.5*(V0b.array().log().sum() + yP0y/vv + (n-p)*log(vv));
    }

    void df_ddf(scalar_t h2, scalar_t & df, scalar_t & ddf) {
      update(h2);
      VECTOR DeltaP0y = Deltab.cwiseProduct( P0y );
      scalar_t y_PDP_y = P0y.dot( DeltaP0y );
  
      VECTOR P0DeltaP0y  = V0bi.asDiagonal() * DeltaP0y;
      scalar_t y_PDPDP_y = DeltaP0y.dot( P0DeltaP0y );

      scalar_t trace_PD = Deltab.cwiseProduct(V0bi).sum();
      scalar_t trace_PDPD = Deltab.cwiseProduct(V0bi.cwiseProduct(Deltab.cwiseProduct(V0bi))).sum();
  
      df = -0.5*(trace_PD - (n-p)*y_PDP_y/yP0y);
      ddf = -0.5*(-trace_PDPD + (n-p)*( 2*y_PDPDP_y/yP0y - (y_PDP_y*y_PDP_y)/(yP0y*yP0y)));
    }

    void blup(scalar_t h2, VECTOR & beta, VECTOR & omega, bool updateh2) {
      if(updateh2) update(h2);
      VECTOR Sigmab = Sigma.bottomRows(n-p);
      omega = h2 * Sigmab.asDiagonal() * P0y;
      beta = Y.topRows(p);
    }

};
#endif
