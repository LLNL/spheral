//---------------------------------Spheral++----------------------------------//
// CubicHermiteInterpolator
//
// An (optionally monotonic) form of cubic Hermite interpolation.
//
// Created by JMO, Fri Apr  1 14:22:04 PDT 2022
//----------------------------------------------------------------------------//
#ifndef __Spheral_CubicHermiteInterpolator__
#define __Spheral_CubicHermiteInterpolator__

#include <cstddef>
#include <vector>

namespace Spheral {

class CubicHermiteInterpolator {
public:
  //--------------------------- Public Interface ---------------------------//
  // Constructors, destructors
  template<typename Func>
  CubicHermiteInterpolator(const double xmin,
                           const double xmax,
                           const size_t n,
                           const Func& F);
  template<typename Func, typename GradFunc>
  CubicHermiteInterpolator(const double xmin,
                           const double xmax,
                           const size_t n,
                           const Func& F,
                           const GradFunc& Fgrad);
  CubicHermiteInterpolator(const double xmin,
                           const double xmax,
                           const std::vector<double>& values);
  CubicHermiteInterpolator();
  ~CubicHermiteInterpolator();

  // Copy and assignment
  CubicHermiteInterpolator(const CubicHermiteInterpolator& rhs);
  CubicHermiteInterpolator& operator=(const CubicHermiteInterpolator& rhs);

  // (Re)initialize after construction, same options as construction
  template<typename Func>
  void initialize(const double xmin,
                  const double xmax,
                  const size_t n,
                  const Func& F);
  template<typename Func, typename GradFunc>
  void initialize(const double xmin,
                  const double xmax,
                  const size_t n,
                  const Func& F,
                  const GradFunc& Fgrad);
  void initialize(const double xmin,
                  const double xmax,
                  const std::vector<double>& yvals);

  // Force interpolation to be monotonic (may introduce structure between tabulated points)
  void makeMonotonic();

  // Comparisons
  bool operator==(const CubicHermiteInterpolator& rhs) const;

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

  // Compute the Hermite basis functions
  double h00(const double x) const;
  double h10(const double x) const;
  double h01(const double x) const;
  double h11(const double x) const;

  // Allow read access the internal data representation
  size_t size() const;                        // The number of tabulated values
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

  // Initialize the gradient at the interpolation points based on the tabulated
  // interpolation values
  void initializeGradientKnots();
};

}

#include "CubicHermiteInterpolatorInline.hh"

#endif
