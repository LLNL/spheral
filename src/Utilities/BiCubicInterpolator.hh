//---------------------------------Spheral++----------------------------------//
// BiCubicInterpolator
//
// Encapsulates the algorithm and data for cubic interpolation in 2D
// Assumes the results is interpolated as
//   <F(x,y)> = [1 x x^2 x^3][c00 c01 c02 c03][1]
//                           [c10 c11 c12 c13][y]
//                           [c20 c21 c22 c23][y^2]
//                           [c30 c31 c32 c33][y^3]
// Assumes we provide functors for F and gradF for fitting, where
//      F(x,y) -> double
//  gradF(x,y) -> Dim<2>::SymTensor
//
// See https://en.wikipedia.org/wiki/Bicubic_interpolation
//
// Created by JMO, Thu Mar 10 16:24:45 PST 2022
//----------------------------------------------------------------------------//
#ifndef __Spheral_BiCubicInterpolator__
#define __Spheral_BiCubicInterpolator__

#include "Utilities/XYInterpolator.hh"

namespace Spheral {

class BiCubicInterpolator: public XYInterpolator {
public:
  //--------------------------- Public Interface ---------------------------//
  // Constructors, destructors
  template<typename Func, typename GradFunc>
  BiCubicInterpolator(const double xmin,
                      const double xmax,
                      const double ymin,
                      const double ymax,
                      const size_t nx,
                      const size_t ny,
                      const Func& F,
                      const GradFunc& gradF,
                      const bool xlog = false,
                      const bool ylog = false);
  BiCubicInterpolator();
  virtual ~BiCubicInterpolator();

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

#include "BiCubicInterpolatorInline.hh"

#else

// Forward declaration
namespace Spheral {
  class BiCubicInterpolator;
}

#endif
