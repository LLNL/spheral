//---------------------------------Spheral++----------------------------------//
// BiLinearInterpolator
//
// Encapsulates the algorithm and data for bi-linear interpolation in 2D
// Assumes the results is interpolated as
//   <F(x,y)> = c0 + c1*x + c2*y + c3*x*y
//
// Created by JMO, Tue Mar  1 10:59:20 PST 2022
//----------------------------------------------------------------------------//
#ifndef __Spheral_BiLinearInterpolator__
#define __Spheral_BiLinearInterpolator__

#include "Utilities/XYInterpolator.hh"

namespace Spheral {

class BiLinearInterpolator: public XYInterpolator {
public:
  //--------------------------- Public Interface ---------------------------//
  // Constructors, destructors
  template<typename Func>
  BiLinearInterpolator(const double xmin,
                       const double xmax,
                       const double ymin,
                       const double ymax,
                       const size_t nx,
                       const size_t ny,
                       const Func& F,
                       const bool xlog,
                       const bool ylog);
  BiLinearInterpolator();
  virtual ~BiLinearInterpolator();

  // Interpolate for the F(x,y) value
  double operator()(const double x, const double y) const;

  // Interpolated gradient values.
  double prime_x(const double x, const double y) const;
  double prime_y(const double x, const double y) const;

  // Return the lower bound index in the table of coefficients for the given position
  size_t lowerBound(const double x, const double y) const;
};

}

#include "BiLinearInterpolatorInline.hh"

#endif
