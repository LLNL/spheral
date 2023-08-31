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

#include <cstddef>
#include <vector>

#include "Field/SphArray.hh"

namespace Spheral {

class QuadraticInterpolator {
  using DataType = double;
  using ContainerType = ManagedVector<DataType>;
public:
  //--------------------------- Public Interface ---------------------------//
  // Constructors, destructors
  QuadraticInterpolator();

  template<typename Func>
  QuadraticInterpolator(const double xmin,
                        const double xmax,
                        const size_t n,
                        const Func& F);

  // Alternatively initialize from tabulated values
  QuadraticInterpolator(const double xmin, const double xmax,
                  const std::vector<double>& yvals);

  RAJA_HOST_DEVICE
  QuadraticInterpolator(const QuadraticInterpolator& qInt);

  // Comparisons
  RAJA_HOST_DEVICE
  bool operator==(const QuadraticInterpolator& rhs) const;

  // Interpolate for the y value
  RAJA_HOST_DEVICE
  double operator()(const double x) const;
  RAJA_HOST_DEVICE
  double prime(const double x) const;    // First derivative
  RAJA_HOST_DEVICE
  double prime2(const double x) const;   // Second derivative

  // Same as above, but use a pre-computed table position (from lowerBound)
  RAJA_HOST_DEVICE
  double operator()(const double x, const size_t i0) const;
  RAJA_HOST_DEVICE
  double prime(const double x, const size_t i0) const;    // First derivative
  RAJA_HOST_DEVICE
  double prime2(const double x, const size_t i0) const;   // Second derivative

  // Return the lower bound index in the table for the given x coordinate
  RAJA_HOST_DEVICE
  size_t lowerBound(const double x) const;

  // Allow read access the internal data representation
  RAJA_HOST_DEVICE
  size_t size() const;                        // The size of the tabulated coefficient arrays
  RAJA_HOST_DEVICE
  double xmin() const;                        // Minimum x coordinate for table              
  RAJA_HOST_DEVICE
  double xmax() const;                        // Maximum x coordinate for table              
  RAJA_HOST_DEVICE
  double xstep() const;                       // delta x between tabulated values            

  RAJA_HOST_DEVICE
  const ContainerType& coeffs() const;  // the fitting coefficients
  
private:
  //--------------------------- Private Interface --------------------------//
  // Member data
  size_t mN1;
  double mXmin, mXmax, mXstep;
  ContainerType mcoeffs;
};

}

#include "QuadraticInterpolatorInline.hh"

#else

// Forward declaration
namespace Spheral {
  class QuadraticInterpolator;
}

#endif
