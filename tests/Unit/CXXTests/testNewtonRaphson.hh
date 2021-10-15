//------------------------------------------------------------------------------
// testNewtonRaphson
// Simple testing of the Newton-Raphson root finder.
//
// Created by JMO, Fri Sep 24 09:37:52 PDT 2004
//------------------------------------------------------------------------------
#ifndef __Spheral_testNewtonRaphson_hh__
#define __Spheral_testNewtonRaphson_hh__

#include "Utilities/newtonRaphson.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Define a functor class which we'll pass to the Newton-Raphson algorithm
// as our function object.  We'll make it a third order polynomial.
//------------------------------------------------------------------------------
class TestFunction {
 public:
  TestFunction(const double x0,
               const double x1,
               const double x2):
    mx0(x0),
    mx1(x1),
    mx2(x2) {}
  ~TestFunction() {}
  std::pair<double, double> operator()(const double x) const {
    return std::pair<double, double>((x - mx0)*(x - mx1)*(x - mx2),
                                     3.0*x*x - 2.0*(mx0 + mx1 + mx2)*x + (mx0*mx1 + mx0*mx2 + mx1*mx2));
  }
 private:
  double mx0, mx1, mx2;
};

//------------------------------------------------------------------------------
// A function we can use to call the Newton-Raphson method with instances of the
// TestFunction class.
//------------------------------------------------------------------------------
double
testNewtonRaphsonRoot(const TestFunction& func,
                      const double x1,
                      const double x2) {
  return newtonRaphson(func, x1, x2);
}

}

#endif
