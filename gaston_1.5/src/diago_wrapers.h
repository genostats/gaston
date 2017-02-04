#include <Rcpp.h>
#include <RcppEigen.h>
#include <iostream>
using namespace Eigen;
using namespace Rcpp;
typedef Map<MatrixXd> Map_MatrixXd;

// structure contenant tous les paramètres pour diago_likelihood
struct par_li {
  int p;
  const MatrixXd * y;
  MatrixXd * x;
  const Map_MatrixXd * sigma;
  VectorXd * P0y;
  double * v;
  MatrixXd * XViXi;
  double likelihood;
};

//wraper pour diago likelihood
double wrap_li(double h2, void * par);
double wrap_li_nc(double h2, void * par);

//wraper pour diago full likelihood
double wrap_full_li(double h2, void * par);
double wrap_full_li_nc(double h2, void * par);

