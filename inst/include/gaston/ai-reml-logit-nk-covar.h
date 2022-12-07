// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <math.h>
#include <iostream>
#include "logit_model.h"
#include "matrix-varia.h"
#include "any_nan.h"

#ifndef GASTONAIREMLn_logit
#define GASTONAIREMLn_logit

//#define ANY(_X_) (std::any_of(_X_.begin(), _X_.end(), [](bool x) {return x;})) 

using namespace Rcpp;
using namespace Eigen;

typedef Map<MatrixXd> Map_MatrixXd;


// AI REML avec n matrices de Kinship
template<typename T1, typename T2, typename A, typename T3, typename T4>

void AIREMLn_logit(const Eigen::MatrixBase<T1> & y, const Eigen::MatrixBase<T4> & x, const std::vector<T2,A> & K, bool constraint, 
               const Eigen::MatrixBase<T3> & min_tau, int max_iter, double eps, bool verbose, VectorXd & tau,
			   int & niter, MatrixXd & P, VectorXd & omega, VectorXd & beta, MatrixXd & XViX_i,
			   bool & start_tau, bool & start_beta, bool EM) {

  int n(y.rows()), p(x.cols()), s(K.size()), i(0);
  int j;
  
  MatrixXd W(n,n), V(n,n), Vi(n,n);
  MatrixXd XViX(p,p);
  MatrixXd ViX(n,p);
  VectorXd z(n), pi(n), Pz(n), PPz(n);

  std::vector<VectorXd> KPz, PKPz;
  for(int i = 0; i < s; i++) {
    KPz.push_back(VectorXd(n));
    PKPz.push_back(VectorXd(n));
  }

  VectorXd tau0(s), gr(s), dif(s+p), beta0(p);
  MatrixXd AI(s, s), pi_AI(s,s);
  double log_detV, detV, d, ld, log_d, d1, log_d1;

  // X'X
  MatrixXd xtx( MatrixXd(p,p).setZero().selfadjointView<Lower>().rankUpdate( x.transpose() ));
  MatrixXd xtxi(p,p); // et son inverse
  double det_xtx, ldet_xtx;
  MatrixXd xtx0(xtx);
  sym_inverse(xtx0, xtxi, ldet_xtx, det_xtx, 1e-5); // détruit xtx0
  
  // initialisation beta
  if (!start_beta) logistic_model2<double>(y, x, beta, XViX_i, 1e-3);
  
  if(verbose) Rcout << "[Initialization] beta = " << beta.transpose() << "\n";    
   
  // initialisation pseudo réponse
  W.setZero();
  for (int j = 0; j < n; j++) {
    pi(j) = 1/( 1 + exp( - x.row(j).dot(beta) ) );
    W(j,j) = 1/( pi(j)*(1-pi(j)) );
    z(j) = x.row(j).dot(beta) + (y(j)-pi(j))/(pi(j)*(1-pi(j))); }
  if (!start_tau) {
    for (int j = 0; j < s; j++) 
      tau(j) = (z.dot(z)/(n-1)- z.sum()*z.sum()/(n-1)/n)/s; 
  }
  if(verbose) Rcout << "[Initialization] tau = " << tau.transpose() << "\n";    

  // first update tau
  V.noalias() = W;
  for (int j = 0; j < s; j++) V.noalias() += tau(j)*K[j];

  sym_inverse(V,Vi,log_detV,detV,1e-7);
  ViX.noalias() = Vi * x;
  XViX.noalias() = x.transpose() * ViX;
  sym_inverse(XViX, XViX_i, log_d, d, 1e-5);
  P.noalias() = Vi - ViX * XViX_i * ViX.transpose();
  Pz.noalias()   =  P.selfadjointView<Lower>() * z;
  
  for(int j = 0; j < s; j++) {
    KPz[j].noalias() = K[j] * Pz;
    gr(j) = -0.5*(trace_of_product(K[j], P) - Pz.dot(KPz[j]));
    tau(j) += 2*tau(j)*tau(j)*gr(j)/n;
  }
  if(verbose) Rcout << "[Iteration " << i+1 << "] tau = " << tau.transpose() << "\n";    
	
  // update omega &  beta
  omega.setZero();
  for(int j = 0; j < s; j++) omega.noalias() += tau(j)*KPz[j];
  beta0 = beta;
  beta.noalias() = x.transpose() * (z - omega - W*Pz);
  beta = xtxi * beta;
  if(verbose) Rcout << "[Iteration " << i+1 << "] beta = " << beta.transpose() << "\n";    

  for(i = 1; i < max_iter; i++) {  
    // update pseudo reponse
    W.setZero();
    for(int j = 0; j < n; j++) {
      pi(j) = 1/( 1 + exp( - x.row(j).dot(beta) - omega(j) ) );
      W(j,j) = 1/( pi(j)*(1-pi(j)) );
      z(j) = x.row(j).dot(beta) + omega(j) + (y(j)-pi(j))/(pi(j)*(1-pi(j))); }
 	
    V.noalias() = W;	
    for(int j = 0; j < s; j++)
      V.noalias() += tau(j)*K[j];

   // Calcul de Vi = inverse(V)
    sym_inverse(V,Vi,log_detV,detV,1e-7);

    // if(verbose) Rcout << "[Iteration " << i+1 << "] det V = " << detV << "\n";


    // Calcul de P
    ViX.noalias() = Vi * x;
    XViX.noalias() = x.transpose() * ViX;
    sym_inverse(XViX, XViX_i, log_d1, d1, 1e-5);
    P.noalias() = Vi - ViX * XViX_i * ViX.transpose();

    // gradient
    Pz.noalias()   =  P.selfadjointView<Lower>() * z;
    for(int j = 0; j < s; j++) {
      KPz[j].noalias() = K[j] * Pz;
      if(!EM) PKPz[j].noalias() = P.selfadjointView<Lower>() * KPz[j];
      gr(j) = -0.5*(trace_of_product(K[j], P) - Pz.dot(KPz[j]));
    }
    //  if(verbose) Rcout << "[Iteration " << i+1 << "] gradient = " << gr.transpose() << "\n";


    // UPDATE tau
    tau0 = tau;
    if(!EM) {
      // updating tau with AIREML
      // Compute Average Information
      AI.setZero();
      for(int j = 0; j < s; j++)
        AI(j,j) = 0.5*PKPz[j].dot(KPz[j]);
      for(int j1 = 1; j1 < s; j1++) {
        for(int j2 = 0; j2 < j1; j2++) {
          AI(j1,j2) = AI(j2,j1) = 0.5*PKPz[j1].dot(KPz[j2]);
        }
      }
      // update tau
      sym_inverse(AI,pi_AI,d,ld,1e-5);
      tau += pi_AI * gr;
    }
    // did it fail?
    bool failed = any_nan(tau);
    if(failed) {
      if(verbose) Rcout << "[Iteration " << i+1 << "] AIREML step failed, falling back to an EM step\n";  
      tau = tau0;
    }
    if(EM || failed) { 
      // updating tau with EM
      for (j = 0; j < s; j++)
        tau(j) += 2*tau(j)*tau(j)*gr(j)/n;
    }

    if(constraint) {
      for(int j = 0; j < s; j++) {
        if(tau(j) < min_tau(j)) {
          tau(j) = min_tau(j);
          if(verbose) Rcout << "[Iteration " << i+1 << "] Constraining tau[" << j << "]\n";
        }
      }
    }
    if(verbose) Rcout << "[Iteration " << i+1 << "] tau = " << tau.transpose() << "\n";  

    // update BLUP
    // BLUP pour omegaz et omega
    omega.setZero();
    for(int j = 0; j < s; j++)
    omega.noalias() += tau(j)*KPz[j];
	  
    // BLUP pour beta
    beta0=beta;
    beta = x.transpose() * (z - omega - W*Pz);
    beta = xtxi * beta;	
    if(verbose) Rcout << "[Iteration " << i+1 << "] beta = " << beta.transpose() << "\n";    

    // Parametre de convergence
    for(int j = 0; j < p; j++)
      dif(j) = fabs(beta0(j)-beta(j))/(fabs(beta0(j))+fabs(beta(j)));
    for(int j = 0; j < s; j++)
      dif(j+p) = fabs(tau0(j)-tau(j))/(fabs(tau0(j))+fabs(tau(j)));
    if(verbose) Rcout << "[Iteration " << i+1 << "] dif = " << dif.transpose() << "\n";    

    if(2*dif.maxCoeff() < eps)
      break; 

    R_CheckUserInterrupt();
  }
  niter = i+1;
  
  V.noalias() = W;
  for(j = 0; j < s; j++)
    V.noalias() += tau(j)*K[j];
  sym_inverse(V,V,log_detV,detV,1e-7);
  ViX.noalias() = Vi * x;
  XViX.noalias() = x.transpose() * ViX;
  sym_inverse(XViX, XViX_i, log_d1, d1, 1e-5);
  P.noalias() = Vi - ViX * XViX_i * ViX.transpose();
  Pz.noalias()   =  P.selfadjointView<Lower>() * z;
  
  omega.setZero();
  for (j=0; j<s; j++) {
    KPz[j].noalias()  = K[j] * Pz;
    omega.noalias() += tau(j)*KPz[j]; 
  }

  beta = x.transpose() * (z - omega - W*Pz);
  beta = xtxi * beta;
}
#endif
