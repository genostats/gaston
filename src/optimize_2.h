#include <iostream>
class fun {
    double scale;
  public:
    fun() : scale(1.0) {};
    virtual double f(double x) { 
      std::cout << "oups";
      return 0; 
    }
    double Brent_fmin(double ax, double bx, double tol);
    double Brent_fmax(double ax, double bx, double tol) {
      scale *= -1;
      double x = Brent_fmin(ax, bx, tol);
      scale *= -1;
      return x;
    }
};

// double Brent_fmin(double ax, double bx, fun * A, double tol);


