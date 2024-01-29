//---------------------------------Spheral++----------------------------------//
// BiQuadraticInterpolator
//
// Encapsulates the algorithm and data for parabolic interpolation in 2D
// Assumes the results is interpolated as
//   <F(x,y)> = c0 + c1*x + c2*y + c3*x*y + c4*x^2 + c5*x^2*y + c6*y^2 + c7*x*y^2 + c8*x^2*y^2
//
// Created by JMO, Thu Dec 10 14:48:01 PST 2020
//----------------------------------------------------------------------------//
#ifndef __Spheral_BiQuadraticInterpolator__
#define __Spheral_BiQuadraticInterpolator__

#include "Utilities/XYInterpolator.hh"

namespace Spheral {

class BiQuadraticInterpolator: public XYInterpolator {
public:
  //--------------------------- Public Interface ---------------------------//
  // Constructors, destructors
  template<typename Func>
  BiQuadraticInterpolator(const double xmin,
                          const double xmax,
                          const double ymin,
                          const double ymax,
                          const size_t nx,
                          const size_t ny,
                          const Func& F,
                          const bool xlog = false,
                          const bool ylog = false);
  BiQuadraticInterpolator();
  virtual ~BiQuadraticInterpolator();

  // Interpolate for the F(x,y) value
  double operator()(const double x, const double y) const;

  // Interpolated gradient values.
  double prime_x(const double x, const double y) const;
  double prime_y(const double x, const double y) const;
  double prime2_xx(const double x, const double y) const;
  double prime2_xy(const double x, const double y) const;
  double prime2_yy(const double x, const double y) const;
};

}

#include "BiQuadraticInterpolatorInline.hh"

#endif
