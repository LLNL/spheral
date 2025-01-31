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

namespace Spheral {

class QuadraticInterpolator {
public:
  //--------------------------- Public Interface ---------------------------//
  // Constructors, destructors
  template<typename Func>
  QuadraticInterpolator(double xmin, double xmax, size_t n, const Func& F);
  QuadraticInterpolator(double xmin, double xmax, const std::vector<double>& yvals);
  QuadraticInterpolator();
  ~QuadraticInterpolator();

  // Initialize after construction, either with a function or tabulated values
  template<typename Func>
  void initialize(double xmin, double xmax, size_t n, const Func& f);
  void initialize(double xmin, double xmax, const std::vector<double>& yvals);

  // Comparisons
  bool operator==(const QuadraticInterpolator& rhs) const;

  // Interpolate for the y value
  double operator()(const double x) const;
  double prime(const double x) const;    // First derivative
  double prime2(const double x) const;   // Second derivative

  // Same as above, but use a pre-computed table position (from lowerBound)
  double operator()(const double x, const size_t i0) const;
  double prime(const double x, const size_t i0) const;    // First derivative
  double prime2(const double x, const size_t i0) const;   // Second derivative

  // Return the lower bound index in the table for the given x coordinate
  size_t lowerBound(const double x) const;

  // Allow read access the internal data representation
  size_t size() const;                        // The size of the tabulated coefficient arrays
  double xmin() const;                        // Minimum x coordinate for table              
  double xmax() const;                        // Maximum x coordinate for table              
  double xstep() const;                       // delta x between tabulated values            
  const std::vector<double>& coeffs() const;  // the fitting coefficients
  
private:
  //--------------------------- Private Interface --------------------------//
  // Member data
  size_t mN1;
  double mXmin, mXmax, mXstep;
  std::vector<double> mcoeffs;
};

}

#include "QuadraticInterpolatorInline.hh"

#endif
