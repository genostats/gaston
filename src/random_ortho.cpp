#include <Rcpp.h>
#include <RcppEigen.h>
using namespace Rcpp;
using namespace Eigen;

// [[Rcpp::export]]
NumericMatrix random_ortho(int n) {
  NumericMatrix R(n,n);
  Map<MatrixXd> r(as<Map<MatrixXd> >(R));
  // matrice 2x2
  double theta = 2*M_PI*R::runif(0,1);
  r(0,0) = cos(theta);
  r(1,1) = cos(theta);
  r(0,1) = sin(theta);
  r(1,0) = -sin(theta);
  for(int i = 2; i < n; i++) {
    // nouvelle colonne aléatoire
    Block<Map<MatrixXd> > V = r.block(0,i,i+1,1);
    for(int j = 0; j <= i; j++) V(j,0) = R::rnorm(0,1);
    V.col(0).normalize();
    // Householder pour transformer les cols précédentes
    VectorXd W = V.col(0);
    W(i) -= 1;
    Block<Map<MatrixXd> > a = r.topLeftCorner(i+1,i);
    // le .eval() permet d'assurer que le membre de droite du produit
    // est évalué dans un objet temporaire
    // cela permet le .noalias() qui impose que la soustraction se fait 'in place'
    a.noalias() -= 2*W*(((W.transpose()/W.squaredNorm())*a).eval());
  }
  return R;
}

RcppExport SEXP gg_random_ortho(SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(random_ortho(n));
    return rcpp_result_gen;
END_RCPP
}

