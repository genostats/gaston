// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <math.h>
#include <iostream>
#include "any_nan.h"
#ifndef GASTONAIREMLn_logit_nofix
#define GASTONAIREMLn_logit_nofix

using namespace Rcpp;
using namespace Eigen;

typedef Map<MatrixXd> Map_MatrixXd;


// AI REML avec n matrices de Kinship
template<typename T1, typename T2, typename A, typename T3>
void AIREMLn_logit_nofix(const Eigen::MatrixBase<T1> & y, const std::vector<T2,A> & K, bool constraint, 
                         const Eigen::MatrixBase<T3> & min_tau, int max_iter, 
                         double eps, bool verbose, VectorXd & tau, int & niter, MatrixXd & P, VectorXd & omega, 
                         bool start_tau, bool EM) {

  int n(y.rows()), s(K.size()), i(0);
  int j;  
  MatrixXd W(n,n), V(n,n);
  VectorXd pi(n), z(n), Pz(n), PPz(n), dif(s);
  
  std::vector<VectorXd> KPz, PKPz;
  for(int i = 0; i < s; i++) {
    KPz.push_back(VectorXd(n));
    PKPz.push_back(VectorXd(n));
  }

  VectorXd tau0(s), gr(s);
  MatrixXd AI(s, s), pi_AI(s,s);
  double log_detV, detV, d, ld;

  // initialisation pseudo réponse
  W.setZero();
  for(int j = 0; j < n; j++) {
    W(j,j) = 4;
    z(j)= (y(j)-1/2)*4; }
	
  if(!start_tau) {
    for (j=0; j<s; j++) tau(j)= ( z.dot(z)/(n-1)- z.sum()*z.sum()/(n-1)/n )/s;
  }  
  if(verbose) Rcout << "[Initialization] tau = " << tau.transpose() << "\n";    

  // Variance matrix
  V.noalias() = W;
  for (j = 0; j < s; j++)
    V.noalias() += tau(j)*K[j];
  sym_inverse(V,P,log_detV,detV,1e-7);
 
  // first update
  Pz.noalias()   =  P.selfadjointView<Lower>() * z;
  for (j = 0; j < s; j++) {
    KPz[j].noalias()  = K[j] * Pz;
    gr(j) = -0.5*(trace_of_product(K[j],P) - Pz.dot(KPz[j])); 
	tau(j) += 2*tau(j)*tau(j)*gr(j)/n;
  }
  if(verbose) Rcout << "[Iteration " << i+1 << "] gr = " << gr.transpose() << "\n";   
  if(verbose) Rcout << "[Iteration " << i+1 << "] tau = " << tau.transpose() << "\n";    
	
  // update omega
  omega.setZero();
  for (int j = 0; j < s; j++)
      omega.noalias() += tau(j)*KPz[j];

  for(i = 1; i < max_iter; i++) {  
    // update pseudo reponse
    W.setZero();
    for(int j = 0; j < n; j++) {
      pi(j) = 1/( 1 + exp( - omega(j) ) );
      W(j,j) = 1/( pi(j)*(1-pi(j)) );
      z(j) = omega(j) + (y(j)-pi(j))/(pi(j)*(1-pi(j))); }
    V.noalias() = W;	  
    for(int j = 0; j < s; j++)
      V.noalias() += tau(j)*K[j];
	  
    // calcul de P = inverse(V)
    sym_inverse(V,P,log_detV, detV,1e-7);
    Pz.noalias()   =  P.selfadjointView<Lower>() * z;
	
    // gradient
    for(int j = 0; j < s; j++) {
      KPz[j].noalias() = K[j] * Pz;
      if(!EM) PKPz[j].noalias() = P.selfadjointView<Lower>() * KPz[j];
      gr(j) = -0.5*(trace_of_product(K[j], P) - Pz.dot(KPz[j]));
    }

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
	
    for(int j = 0; j < s; j++)
      dif(j) = fabs( tau0(j)-tau(j) )/( fabs(tau0(j))+ fabs(tau(j)) );
	  	
    if(2*dif.maxCoeff() < eps) break; 

    R_CheckUserInterrupt();
	
    // calcul de omega (breeding value)
    omega.setZero();
    for(int j = 0; j < s; j++)
      omega.noalias() += tau(j)*KPz[j];
  }
	
  V.noalias() = W;
  for(int j = 0; j < s; j++)
    V.noalias() += tau(j)*K[j];
  sym_inverse(V,P,log_detV,detV,1e-7);
  Pz.noalias()   =  P.selfadjointView<Lower>() * z;
  
  omega.setZero();
  for (j=0; j<s; j++) {
  KPz[j].noalias()  = K[j] * Pz;
  omega.noalias() += tau(j)*KPz[j]; }
  
  niter = i+1;
}
#endif
