#include <RcppEigen.h>
#include <math.h>
#include <iostream>
#include "matrix-varia.h"
#include "logit.h"
 
using namespace Rcpp;
using namespace Eigen;

// float [pas besoin de template, on n'aura pas de MapMatrix]
void AIREML1_logit_f(const MatrixXf & y, const MatrixXf & x, const MatrixXf & K, 
                     bool constraint, float min_tau, int max_iter, float eps, bool verbose, float & tau, int & niter, MatrixXf & P,
                     VectorXf & omega, VectorXf & beta, MatrixXf & XViX_i, bool start_tau, bool start_beta) {

  int n(y.rows()), p(x.cols()), i(0);
  MatrixXf V(n,n), Vi(n,n);
  VectorXf W(n);
  MatrixXf XViX(p,p), ViX(n,p);
  VectorXf dif(p+1), beta0(p);
  VectorXf pi(n), z(n), Pz(n), KPz(n), PKPz(n);
  float tau0, log_detV, detV, d, log_d, AI, gr;
  // X'X
  MatrixXf xtx( MatrixXf(p,p).setZero().selfadjointView<Lower>().rankUpdate( x.transpose() ));
  MatrixXf xtxi(p,p); // et son inverse
  float det_xtx, ldet_xtx;
  MatrixXf xtx0(xtx);
  sym_inverse(xtx0, xtxi, ldet_xtx, det_xtx, 1e-5); // détruit xtx0
  
  // initialisation beta
  if(!start_beta) logistic_model_f(y, x, 1e-3, beta, XViX_i); 
  
  if(verbose) Rcout << "[Initialization] beta = " << beta.transpose() << "\n";    
   
  // initialisation pseudo réponse
  // W.setZero();
  for(int j = 0; j < n; j++) {
    pi(j) = 1/( 1 + exp( - x.row(j).dot(beta) ) );
    W(j) = 1/( pi(j)*(1-pi(j)) );
    z(j) = x.row(j).dot(beta) + (y(j)-pi(j))/(pi(j)*(1-pi(j)));
  }
  if (!start_tau) tau= z.dot(z)/(n-1)- z.sum()*z.sum()/(n-1)/n;
  if(verbose) Rcout << "[Initialization] tau = " << tau << "\n";    

  // first update tau
  // V.noalias() = W.asDiagonal() + tau*K;
  V.noalias() = tau*K; V.diagonal().noalias() += W;
  sym_inverse(V,Vi,log_detV,detV,1e-7);
  ViX.noalias() = Vi * x;
  XViX.noalias() = x.transpose() * ViX;
  sym_inverse(XViX, XViX_i, log_d, d, 1e-5);
  P.noalias() = Vi - ViX * XViX_i * ViX.transpose();
  Pz.noalias()   =  P.selfadjointView<Lower>() * z;
  KPz.noalias()  = K * Pz; // le .selfadjointView ne compile pas avec le template !!  
  
  gr = -0.5*(trace_of_product(K,P) - Pz.dot(KPz));   
  tau += 2*tau*tau*gr/n;
  if(verbose) Rcout << "[Iteration " << i+1 << "] tau = " << tau << "\n";    
	
  // update omega &  beta
  omega.noalias() = tau*KPz;
  beta0 = beta;
  beta.noalias() = x.transpose() * (z - omega - W.asDiagonal()*Pz);
  beta = xtxi * beta;
  if(verbose) Rcout << "[Iteration " << i+1 << "] beta = " << beta.transpose() << "\n";    
	  
  // Itération
  for(i = 1; i < max_iter; i++) {
    // update pseudo reponse
    for(int j = 0; j < n; j++) {
      pi(j) = 1/( 1 + exp( - x.row(j).dot(beta) - omega(j) ) );
      W(j) = 1/( pi(j)*(1-pi(j)) );
      z(j) = x.row(j).dot(beta) + omega(j) + (y(j)-pi(j))/(pi(j)*(1-pi(j)));
    }
	  
    // update P
    V.noalias() = tau*K; V.diagonal().noalias() += W;
    sym_inverse(V,Vi,log_detV,detV,1e-7);
    ViX.noalias() = Vi * x;
    XViX.noalias() = x.transpose() * ViX;
    sym_inverse(XViX, XViX_i, log_d, d, 1e-5);
    P.noalias() = Vi - ViX * XViX_i * ViX.transpose();
 
    // gradient
    Pz.noalias()   =  P.selfadjointView<Lower>() * z;
    KPz.noalias()  = K * Pz; // le .selfadjointView ne compile pas avec le template !!
    gr = -0.5*(trace_of_product(K,P) - Pz.dot(KPz));
    if(verbose) Rcout << "[Iteration " << i+1 << "] gr = " << gr << "\n";
	
    //update tau
    // Average Information
    tau0 = tau;
    PKPz.noalias() = P.selfadjointView<Lower>() * KPz;   
    AI = 0.5*PKPz.dot(KPz);
    tau += gr/AI;
    
    if(constraint && tau < min_tau) {
      tau = min_tau;
      if(verbose) Rcout << "[Iteration " << i+1 << "] Constraining tau\n";
    }	
    if(verbose) Rcout << "[Iteration " << i+1 << "] tau = " << tau << "\n";    

    // update omega & beta
    omega.noalias() = tau*KPz;
    beta0 = beta;
    beta.noalias() = x.transpose() * (z - omega - W.asDiagonal()*Pz);
    beta = xtxi * beta;    
    if(verbose) Rcout << "[Iteration " << i+1 << "] beta = " << beta.transpose() << "\n";    
      
    // Parametre de convergence
    for(int j = 0; j < p; j++)
      dif(j) = fabs(beta0(j)-beta(j))/(fabs(beta0(j))+fabs(beta(j)));
    dif(p) = fabs(tau0-tau)/(fabs(tau0)+fabs(tau));
        
    if(2*dif.maxCoeff() < eps) break; 

    checkUserInterrupt();
  }
  niter = i+1;
  
  // update
  V.noalias() = tau*K; V.diagonal().noalias() += W;
  sym_inverse(V,Vi,log_detV,detV,1e-7);
  ViX.noalias() = Vi * x;
  XViX.noalias() = x.transpose() * ViX;
  sym_inverse(XViX, XViX_i, log_d, d, 1e-5);
  P.noalias() = Vi - ViX * XViX_i * ViX.transpose();
  Pz.noalias()   =  P.selfadjointView<Lower>() * z;
  KPz.noalias()  = K * Pz; // le .selfadjointView ne compile pas avec le template !!
  omega.noalias() = tau*KPz;
  beta.noalias() = x.transpose() * (z - omega - W.asDiagonal()*Pz);
  beta = xtxi * beta;    
}
