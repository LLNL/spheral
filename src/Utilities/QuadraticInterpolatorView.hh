//---------------------------------Spheral++----------------------------------//
// QuadraticInterpolatorView
//
// Encapsulates the algorithm and data for parabolic interpolation in 1D
// Assumes the results is interpolated as y_interp = a + b*x + c*x^2
//
// Created by JMO, Fri Dec  4 14:28:08 PST 2020
//----------------------------------------------------------------------------//
#ifndef __Spheral_QuadraticInterpolatorView__
#define __Spheral_QuadraticInterpolatorView__

#include <cstddef>
#include <vector>

#include "Field/SphArray.hh"

namespace Spheral {



class QuadraticInterpolatorView{
public:
  //--------------------------- Public Interface ---------------------------//
  // Constructors, destructors
  SPHERAL_HOST_DEVICE QuadraticInterpolatorView();
  SPHERAL_HOST_DEVICE ~QuadraticInterpolatorView();

  SPHERAL_HOST_DEVICE QuadraticInterpolatorView(QuadraticInterpolatorView const& rhs) = default;

  // Comparisons
  SPHERAL_HOST_DEVICE bool operator==(const QuadraticInterpolatorView& rhs) const;

  // Interpolate for the y value
  SPHERAL_HOST_DEVICE double operator()(const double x) const;
  SPHERAL_HOST_DEVICE double prime(const double x) const;    // First derivative
  SPHERAL_HOST_DEVICE double prime2(const double x) const;   // Second derivative

  // Same as above, but use a pre-computed table position (from lowerBound)
  SPHERAL_HOST_DEVICE double operator()(const double x, const size_t i0) const;
  SPHERAL_HOST_DEVICE double prime(const double x, const size_t i0) const;    // First derivative
  SPHERAL_HOST_DEVICE double prime2(const double x, const size_t i0) const;   // Second derivative

  // Return the lower bound index in the table for the given x coordinate
  SPHERAL_HOST_DEVICE size_t lowerBound(const double x) const;

  // Allow read access the internal data representation
  SPHERAL_HOST_DEVICE size_t size() const;                        // The size of the tabulated coefficient arrays
  SPHERAL_HOST_DEVICE double xmin() const;                        // Minimum x coordinate for table              
  SPHERAL_HOST_DEVICE double xmax() const;                        // Maximum x coordinate for table              
  SPHERAL_HOST_DEVICE double xstep() const;                       // delta x between tabulated values            
                                                                  
  //using CoeffsType = std::vector<double>;
  using CoeffsType = MVSmartRef<double>;
  //const std::vector<double>& coeffs() const;  // the fitting coefficients
  const CoeffsType& coeffs() const;  // the fitting coefficients
  
protected:
  //--------------------------- Private Interface --------------------------//
  // Member data
  size_t mN1;
  double mXmin, mXmax, mXstep;
  CoeffsType mcoeffs;

};

}

#include "QuadraticInterpolatorViewInline.hh"

#endif
