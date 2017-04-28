#include <Rcpp.h>
#include <iostream>
#ifndef GASTON_OPTIM
#define GASTON_OPTIM
template<typename scalar_t>
class fun {
    scalar_t scale;
    scalar_t F(scalar_t x) {
      return scale*f(x);
    }
    void DF_DDF(scalar_t x, scalar_t & df, scalar_t & ddf) {
      df_ddf(x, df, ddf);
      df *= scale; 
      ddf *= scale;
    }

  public:
    fun() : scale(1.0) {};

    virtual scalar_t f(scalar_t x) { 
      return 0; 
    }
    virtual void df_ddf(scalar_t x, scalar_t & df, scalar_t & ddf) {
      df = ddf = 0; 
    }

    scalar_t Brent_fmin(scalar_t ax, scalar_t bx, scalar_t tol);

    scalar_t Brent_fmax(scalar_t ax, scalar_t bx, scalar_t tol) {
      scale = -1;
      scalar_t x = Brent_fmin(ax, bx, tol);
      scale = 1;
      return x;
    }

    void newton_min(scalar_t & x, const scalar_t min_x, const scalar_t max_x, const scalar_t eps, int max_iter, const bool verbose);

    void newton_max(scalar_t & x, const scalar_t min_x, const scalar_t max_x, const scalar_t eps, int max_iter, const bool verbose) {
      scale = -1;
      newton_min(x, min_x, max_x, eps, max_iter, verbose);
      scale = 1;

    }
};


template<typename scalar_t>
void fun<scalar_t>::newton_min(scalar_t & x, const scalar_t min_x, const scalar_t max_x, const scalar_t eps, int max_iter, const bool verbose) {
  int nb_reseeds = 0, i = 0;
  scalar_t df = 1+2*eps;

  bool tried_max = false, tried_min = false;
  if(x == min_x) tried_min = true;
  if(x == max_x) tried_max = true;

  while(std::abs(df) > 2*eps) {
    i++;
    if(i > max_iter) {
      if(verbose) Rcpp::Rcout << "[Iteration " << i << "] Too many iterations, using Brent algorithm" << std::endl;
      x = Brent_fmin(min_x, max_x, 1e-5);
      if(verbose) Rcpp::Rcout << "[Iteration " << i << "] Brent gives " << x << std::endl;
      break;
    }
    scalar_t ddf;
    DF_DDF(x, df, ddf);

    if(verbose) {
      Rcpp::Rcout << "[Iteration " << i << "] ";
      Rcpp::Rcout << "Current point = " << x << " df = " << scale*df << std::endl;
    }

    // si on est au bord de l'intervalle et qu'on ne tend pas à revenir dedans
    if(x == min_x && !std::isnan(df) && df > 0) {
      if(verbose) Rcpp::Rcout << "[Iteration " << i << "] maximum at min = " << x << std::endl;
      break;
    }
    if(x == max_x && !std::isnan(df) && df < 0) {
      if(verbose) Rcpp::Rcout << "[Iteration " << i << "] maximum at max = " << x << std::endl;
      break;
    }

    // si la convexité est mauvaise
    if(ddf < 0) {
      if(verbose) Rcpp::Rcout << "[Iteration " << i << "] likelihood isn't concave" << std::endl;
      scalar_t old_x = x;
      if(df < 0) {
        if(!tried_max) {
          x = max_x;
          tried_max = true;
          if(verbose) Rcpp::Rcout << "[Iteration " << i << "] restarting from " << x << std::endl;
          continue;
        }
        if(verbose) Rcpp::Rcout << "[Iteration " << i << "] Using Brent algorithm" << std::endl;
        x = Brent_fmin(old_x, max_x, 1e-5);
        if(verbose) Rcpp::Rcout << "[Iteration " << i << "] Brent gives " << x << std::endl;
        break;
      }
      if(df > 0) {
        if(!tried_min) {
          x = min_x;
          tried_min = true;
          if(verbose) Rcpp::Rcout << "[Iteration " << i << "] restarting from " << x << std::endl;
          continue;
        }
        if(verbose) Rcpp::Rcout << "[Iteration " << i << "] Using Brent algorithm" << std::endl;
        x = Brent_fmin(min_x, old_x, 1e-5);
        if(verbose) Rcpp::Rcout << "[Iteration " << i << "] Brent gives " << x << std::endl;
        break;
      }
    }
    x -= df/ddf;

    if(std::isnan(x)) {
      if(nb_reseeds++ < 5) {
        x = R::runif(min_x,max_x);
        if(verbose) Rcpp::Rcout << "[Iteration " << i << "] restarting from random value " << x << std::endl;
      } else {
        if(verbose) Rcpp::Rcout << "[Iteration " << i << "] canceling optimization" << std::endl;
        return;
      }
    } else if(x < min_x) {
      x = min_x;
      tried_min = true;
    } else if(x > max_x) {
      x = max_x;
      tried_max = true;
    }
  }
}

// ceci piqué dans le code de R 
// pour utiliser optimize() depuis le code C !
template<typename scalar_t>
scalar_t fun<scalar_t>::Brent_fmin(scalar_t ax, scalar_t bx, scalar_t tol) {
    /*  c is the squared inverse of the golden ratio */
    const scalar_t c = (3. - sqrt(5.)) * .5;

    /* Local variables */
    scalar_t a, b, d, e, p, q, r, u, v, w, x;
    scalar_t t2, fu, fv, fw, fx, xm, eps, tol1, tol3;

/*  eps is approximately the square root of the relative machine precision. */
    eps = DBL_EPSILON;
    tol1 = eps + 1.;/* the smallest 1.000... > 1 */
    eps = sqrt(eps);

    a = ax;
    b = bx;
    v = a + c * (b - a);
    w = v;
    x = v;

    d = 0.;/* -Wall */
    e = 0.;
    fx = F(x);
    fv = fx;
    fw = fx;
    tol3 = tol / 3.;

/*  main loop starts here ----------------------------------- */

    for(;;) {
	xm = (a + b) * .5;
	tol1 = eps * fabs(x) + tol3;
	t2 = tol1 * 2.;

	/* check stopping criterion */

	if (fabs(x - xm) <= t2 - (b - a) * .5) break;
	p = 0.;
	q = 0.;
	r = 0.;
	if (fabs(e) > tol1) { /* fit parabola */
	    r = (x - w) * (fx - fv);
	    q = (x - v) * (fx - fw);
	    p = (x - v) * q - (x - w) * r;
	    q = (q - r) * 2.;
	    if (q > 0.) p = -p; else q = -q;
	    r = e;
	    e = d;
	}

	if (fabs(p) >= fabs(q * .5 * r) || p <= q * (a - x) || p >= q * (b - x)) { /* a golden-section step */
	    if (x < xm) e = b - x; else e = a - x;
	    d = c * e;
	}
	else { /* a parabolic-interpolation step */
	    d = p / q;
	    u = x + d;
	    /* f must not be evaluated too close to ax or bx */
	    if (u - a < t2 || b - u < t2) {
		d = tol1;
		if (x >= xm) d = -d;
	    }
	}

	/* f must not be evaluated too close to x */
	if (fabs(d) >= tol1)
	    u = x + d;
	else if (d > 0.)
	    u = x + tol1;
	else
	    u = x - tol1;

	fu = F(u);
	/*  update  a, b, v, w, and x */

	if (fu <= fx) {
	    if (u < x) b = x; else a = x;
	    v = w;    w = x;   x = u;
	    fv = fw; fw = fx; fx = fu;
	} else {
	    if (u < x) a = u; else b = u;
	    if (fu <= fw || w == x) {
		v = w; fv = fw;
		w = u; fw = fu;
	    } else if (fu <= fv || v == x || v == w) {
		v = u; fv = fu;
	    }
	}
    }
    /* end of main loop */

    return x;
}

#endif
