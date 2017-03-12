// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <iostream>
using namespace Rcpp;
using namespace Eigen;

typedef Map<MatrixXd> Map_MatrixXd;


template<typename T1, typename T2>
inline void chol_inverse(Eigen::MatrixBase<T1> & x, Eigen::MatrixBase<T2> & xi, double & log_det) {
  LDLT<MatrixXd> ldlt(x);

  log_det = ldlt.vectorD().array().log().sum();
  xi.setIdentity();
  ldlt.solveInPlace(xi);
}

// [[Rcpp::export]]
List chol_inverse(NumericMatrix X) {
  Map_MatrixXd x(as<Map<MatrixXd> >(X));
  double log_det;
  NumericMatrix Xi(X.rows(), X.cols());
  Map_MatrixXd xi(as<Map<MatrixXd> >(Xi));
  chol_inverse(x, xi, log_det);
  List L;
  L["inverse"] = Xi;
  L["log_det"] = log_det;
  return L;
}



// cette fonction badibulgue x en l'utilisant pour calculer SD 
// dans le cas où on n'a plus besoin de x une fois qu'on a son inverse ça le fait
// Cette fonction n'utilise que le triangle supérieur de x...
// et ne remplit que le triangle supérieur de y !!
// avec eps = 0 calcule l'inverse
// avec eps petit... (1e-6 ?) pseudo inverse
void blocki(Eigen::MatrixXd & x, int x1, int n, Eigen::MatrixXd & y, int y1, double & log_det, double & det, double eps) {
  if(n == 1) {
    double d = (std::abs(x(x1,x1))<eps)?0:x(x1,x1);
    y(y1,y1) = (d==0)?0:1/d;
    det = d;
    log_det = log(d);
    return;
  }

  int m1 = n/2;
  int m2 = n-m1;

  Block<MatrixXd> A = x.block(x1,x1,m1,m1); 
  //Block<MatrixXd> D = x.block(x1+m1,x1+m1,m2,m2); 
  Block<MatrixXd> B = x.block(x1,x1+m1,m1,m2); 

  Block<MatrixXd> TL = y.block(y1,y1,m1,m1);
  Block<MatrixXd> BL = y.block(y1+m1,y1,m2,m1);
  Block<MatrixXd> TR = y.block(y1,y1+m1,m1,m2);
  Block<MatrixXd> BR = y.block(y1+m1,y1+m1,m2,m2);

  // BR = inverse(D)
  double log_detD, detD;
  blocki(x,x1+m1,m2,y,y1+m1,log_detD, detD, eps);

  // BL = inverse(D)*Bt
  BL.noalias() = BR.selfadjointView<Upper>() * B.transpose();

  // le bloc A est écrasé par SD
  A.triangularView<Upper>() -= B*BL; // on ne calcule que le triangle supérieur car on n'utilise pas l'autre

  // TL = inverse(SD)
  double log_detSD, detSD;
  blocki(x,x1,m1,y,y1,log_detSD, detSD, eps);

  TR.noalias() = TL.selfadjointView<Upper>()*(-BL.transpose());
  BR.triangularView<Upper>() -= BL*TR;  // on sait que cette matrice doit être symmétrique : on ne fait que la moitié des calculs

  log_det = log_detD + log_detSD;
  det = detD * detSD;
}

// la même en float
void blocki(Eigen::MatrixXf & x, int x1, int n, Eigen::MatrixXf & y, int y1, float & log_det, float & det, float eps) {
  if(n == 1) {
    float d = (std::abs(x(x1,x1))<eps)?0:x(x1,x1);
    y(y1,y1) = (d==0)?0:1/d;
    det = d;
    log_det = log(d);
    return;
  }

  int m1 = n/2;
  int m2 = n-m1;

  Block<MatrixXf> A = x.block(x1,x1,m1,m1); 
  //Block<MatrixXf> D = x.block(x1+m1,x1+m1,m2,m2); 
  Block<MatrixXf> B = x.block(x1,x1+m1,m1,m2); 

  Block<MatrixXf> TL = y.block(y1,y1,m1,m1);
  Block<MatrixXf> BL = y.block(y1+m1,y1,m2,m1);
  Block<MatrixXf> TR = y.block(y1,y1+m1,m1,m2);
  Block<MatrixXf> BR = y.block(y1+m1,y1+m1,m2,m2);

  // BR = inverse(D)
  float log_detD, detD;
  blocki(x,x1+m1,m2,y,y1+m1,log_detD, detD, eps);

  // BL = inverse(D)*Bt
  BL.noalias() = BR.selfadjointView<Upper>() * B.transpose();

  // le bloc A est écrasé par SD
  A.triangularView<Upper>() -= B*BL; // on ne calcule que le triangle supérieur car on n'utilise pas l'autre

  // TL = inverse(SD)
  float log_detSD, detSD;
  blocki(x,x1,m1,y,y1,log_detSD, detSD, eps);

  TR.noalias() = TL.selfadjointView<Upper>()*(-BL.transpose());
  BR.triangularView<Upper>() -= BL*TR;  // on sait que cette matrice doit être symmétrique : on ne fait que la moitié des calculs

  log_det = log_detD + log_detSD;
  det = detD * detSD;
}
