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

#include <vector>

namespace Spheral {

class BiLinearInterpolator {
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
                          const Func& F);
  BiLinearInterpolator();
  ~BiLinearInterpolator();

  // Initialize for interpolating the given function
  template<typename Func>
  void initialize(const double xmin,
                  const double xmax,
                  const double ymin,
                  const double ymax,
                  const size_t nx,
                  const size_t ny,
                  const Func& F);

  // Interpolate for the F(x,y) value
  double operator()(const double x, const double y) const;

  // Interpolated gradient values.
  double prime_x(const double x, const double y) const;
  double prime_y(const double x, const double y) const;

  // Return the lower bound index in the table of coefficients for the given position
  size_t lowerBound(const double x, const double y) const;

  // Allow read access the internal data representation
  size_t size() const;                        // The size of the tabulated coefficient arrays
  double xmin() const;                        // Minimum x coordinate for table              
  double xmax() const;                        // Maximum x coordinate for table              
  double ymin() const;                        // Minimum y coordinate for table              
  double ymax() const;                        // Maximum y coordinate for table              
  double xstep() const;                       // x step size
  double ystep() const;                       // y step size
  const std::vector<double>& coeffs() const;  // the fitting coefficients
  
  // Comparison
  bool operator==(const BiLinearInterpolator& rhs) const;

private:
  //--------------------------- Private Interface --------------------------//
  // Member data
  size_t mnx1, mny1;
  double mxmin, mxmax, mymin, mymax, mxstep, mystep;
  std::vector<double> mcoeffs;
};

}

#include "BiLinearInterpolatorInline.hh"

#else

// Forward declaration
namespace Spheral {
  class BiLinearInterpolator;
}

#endif
