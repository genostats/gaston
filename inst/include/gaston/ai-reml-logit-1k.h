// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <iostream>
#ifndef GASTONAIREML1_logit_nofix
#define GASTONAIREML1_logit_nofix
using namespace Rcpp;
using namespace Eigen;

typedef Map<MatrixXd> Map_MatrixXd;

template<typename T1, typename T2>
void AIREML1_logit_nofix(const Eigen::MatrixBase<T1> & y, const Eigen::MatrixBase<T2> & K, bool constraint, double min_tau,
                         int max_iter, double eps, bool verbose, double & tau, int & niter, MatrixXd & P, VectorXd & omega, 
                         bool start_theta, bool EM) {
  int n(y.rows()), i(0);  
  MatrixXd V(n,n), W(n,n);
  VectorXd pi(n), z(n), Pz(n), KPz(n), PKPz(n);
  double tau0, log_detV, detV, AI, gr, dif;
     
  // initialisation pseudo réponse
  W.setZero();
  for(int j = 0; j < n; j++) {
    W(j,j) = 4;
    z(j)= (y(j)-1/2)*4; }
	
  if(!start_theta) {
  tau= z.dot(z)/(n-1)- z.sum()*z.sum()/(n-1)/n;
  }  
  if(verbose) Rcout << "[Initialization] tau = " << tau << "\n";    

  // Variance matrix
  V.noalias() = W + tau*K;
  sym_inverse(V,P,log_detV,detV,1e-7);
 
  // first update
  Pz.noalias()   =  P.selfadjointView<Lower>() * z;
  KPz.noalias()  = K * Pz; // le .selfadjointView ne compile pas avec le template !!
  gr = -0.5*(trace_of_product(K,P) - Pz.dot(KPz)); 
  if(verbose) Rcout << "[Iteration " << i+1 << "] gr = " << gr << "\n";  
  tau += 2*tau*tau*gr/n;
  if(verbose) Rcout << "[Iteration " << i+1 << "] tau = " << tau << "\n";    
	
  // update omega
  omega.noalias() = tau*KPz;
		  
  // Itération
  for(i = 1; i < max_iter; i++) {
    // update pseudo reponse
    for(int j = 0; j < n; j++) {
      pi(j) = 1/( 1 + exp( -omega(j) ) );
      W(j,j) = 1/( pi(j)*(1-pi(j)) );
      z(j) = omega(j) + (y(j)-pi(j))/( pi(j)*(1-pi(j)) ); }
	
    // update P
    V.noalias() = W + tau*K;
    sym_inverse(V,P,log_detV,detV,1e-7);
 
    // gradient
    Pz.noalias()   =  P.selfadjointView<Lower>() * z;
    KPz.noalias()  = K * Pz; // le .selfadjointView ne compile pas avec le template !!
    gr = -0.5*(trace_of_product(K,P) - Pz.dot(KPz));
    if(verbose) Rcout << "[Iteration " << i+1 << "] gr = " << gr << "\n";
	  
    // UPDATE tau
    tau0 = tau;
    if(!EM) {
      // updating tau with AIREML
      // Compute Average Information
      PKPz.noalias() = P.selfadjointView<Lower>() * KPz;   
      AI = 0.5*PKPz.dot(KPz);
      // update tau
      tau += gr/AI;
    }
    // did it fail?
    bool failed = std::isnan(tau);
    if(failed) {
      if(verbose) Rcout << "[Iteration " << i+1 << "] AIREML step failed, falling back to an EM step\n";
      tau = tau0;
    }
    if(EM || failed) {
      // updating tau with EM
      tau += 2*tau*tau*gr/n;
    }
    
    if(constraint && tau < min_tau) {
        tau = min_tau;
        if(verbose) Rcout << "[Iteration " << i+1 << "] Constraining tau\n";
    }	
    if(verbose) Rcout << "[Iteration " << i+1 << "] tau = " << tau << "\n";    

    // update omega
    omega.noalias() = tau*KPz;
          
    // Parametre de convergence
    dif = fabs(tau0-tau)/(fabs(tau0)+fabs(tau));
        
     if(2*dif < eps)
       break;

    checkUserInterrupt();
  }
  niter = i+1;
  
  // update P and omega
  V.noalias() = W + tau*K;
  sym_inverse(V,P,log_detV,detV,1e-7);
  Pz.noalias()   =  P.selfadjointView<Lower>() * z;
  KPz.noalias()  = K * Pz;
  omega.noalias() = tau*KPz;
}
#endif
