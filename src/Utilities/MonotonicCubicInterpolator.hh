//---------------------------------Spheral++----------------------------------//
// MonotonicCubicInterpolator
//
// A monotonic form of cubic Hermite interpolation.
//
// Created by JMO, Fri Apr  1 14:22:04 PDT 2022
//----------------------------------------------------------------------------//
#ifndef __Spheral_MonotonicCubicInterpolator__
#define __Spheral_MonotonicCubicInterpolator__

#include <cstddef>
#include <vector>

namespace Spheral {

class MonotonicCubicInterpolator {
public:
  //--------------------------- Public Interface ---------------------------//
  // Constructors, destructors
  template<typename Func>
  MonotonicCubicInterpolator(const double xmin,
                             const double xmax,
                             const size_t n,
                             const Func& F);
  template<typename Func, typename GradFunc>
  MonotonicCubicInterpolator(const double xmin,
                             const double xmax,
                             const size_t n,
                             const Func& F,
                             const GradFunc& Fgrad);
  MonotonicCubicInterpolator();
  ~MonotonicCubicInterpolator();

  // Alternatively initialize from tabulated values
  void initialize(const double xmin, const double xmax,
                  const std::vector<double>& yvals);

  // Force interpolation to be monotonic (may introduce structure between tabulated points)
  void makeMonotonic();

  // Comparisons
  bool operator==(const MonotonicCubicInterpolator& rhs) const;

  // Interpolate for the y value
  double operator()(const double x) const;

  // Same as above, but use a pre-computed table position (from lowerBound)
  double operator()(const double x, const size_t i0) const;

  // Return the lower bound index in the table for the given x coordinate
  size_t lowerBound(const double x) const;

  // Compute the Hermite basis functions
  double h00(const double x) const;
  double h10(const double x) const;
  double h01(const double x) const;
  double h11(const double x) const;

  // Allow read access the internal data representation
  size_t N() const;                           // The number of tabulated values
  double xmin() const;                        // Minimum x coordinate for table              
  double xmax() const;                        // Maximum x coordinate for table              
  double xstep() const;                       // delta x between tabulated values            
  const std::vector<double>& vals() const;    // the tabulated function values and gradients (sequentially)
  
private:
  //--------------------------- Private Interface --------------------------//
  // Member data
  size_t mN;
  double mXmin, mXmax, mXstep;
  std::vector<double> mVals;
};

}

#include "MonotonicCubicInterpolatorInline.hh"

#else

// Forward declaration
namespace Spheral {
  class MonotonicCubicInterpolator;
}

#endif
