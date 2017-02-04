#include "diago_wrapers.h"
#include "diago.h"
#include "diago_nocovar.h"
#include "diago_full.h"
#include "diago_full_nocovar.h"

//wraper pour diago likelihood
double wrap_li_nc(double h2, void * par) {
  double likelihood = diago_likelihood(h2, ((par_li *) par)->p, *((par_li *) par)->y, *((par_li *) par)->sigma,
                                *((par_li *) par)->P0y, *((par_li *) par)->v);
  ((par_li *) par)->likelihood = likelihood;
  return -likelihood;
}

//wraper pour diago likelihood
double wrap_li(double h2, void * par) {
  double likelihood = diago_likelihood(h2, ((par_li *) par)->p, *((par_li *) par)->y, *((par_li *) par)->x, *((par_li *) par)->sigma,
           *((par_li *) par)->P0y, *((par_li *) par)->v, *((par_li *) par)->XViXi);
  ((par_li *) par)->likelihood = likelihood;
  return -likelihood;
}

//wraper pour diago likelihood
double wrap_full_li(double h2, void * par) {
  double likelihood = diago_full_likelihood(h2, ((par_li *) par)->p, *((par_li *) par)->y, *((par_li *) par)->x, *((par_li *) par)->sigma,
           *((par_li *) par)->P0y, *((par_li *) par)->v, *((par_li *) par)->XViXi);
  ((par_li *) par)->likelihood = likelihood;
  return -likelihood;
}

//wraper pour diago likelihood
double wrap_full_li_nc(double h2, void * par) {
  double likelihood = diago_full_likelihood(h2, ((par_li *) par)->p, *((par_li *) par)->y, *((par_li *) par)->sigma,
           *((par_li *) par)->P0y, *((par_li *) par)->v);
  ((par_li *) par)->likelihood = likelihood;
  return -likelihood;
}

