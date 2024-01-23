//---------------------------------Spheral++----------------------------------//
// QuadraticInterpolator
//
// Encapsulates the algorithm and data for parabolic interpolation in 1D
// Assumes the results is interpolated as y_interp = a + b*x + c*x^2
//
// Created by JMO, Fri Dec  4 14:28:08 PST 2020
//----------------------------------------------------------------------------//
#ifndef __Spheral_QuadraticInterpolator__
#define __Spheral_QuadraticInterpolator__

#include "QuadraticInterpolatorView.hh"

namespace Spheral {

class QuadraticInterpolator : public QuadraticInterpolatorView{
public:
  //--------------------------- Public Interface ---------------------------//
  // Constructors, destructors
  template<typename Func>
  QuadraticInterpolator(const double xmin,
                        const double xmax,
                        const size_t n,
                        const Func& F);
  QuadraticInterpolator() {};
  ~QuadraticInterpolator() {};

  // Alternatively initialize from tabulated values
  void initialize(const double xmin, const double xmax,
                  const std::vector<double>& yvals);

};


}

#endif
