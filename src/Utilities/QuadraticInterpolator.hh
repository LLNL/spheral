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

#include "Utilities/ValueViewInterface.hh"

namespace Spheral {

VVI_IMPL_BEGIN

class QuadraticInterpolator {
public:
  using CoeffsType = vvi::vector<double>;
  //--------------------------- Public Interface ---------------------------//
  // Constructors, destructors
  template<typename Func>
  QuadraticInterpolator(double xmin, double xmax, size_t n, const Func& F);
  QuadraticInterpolator(double xmin, double xmax, const std::vector<double>& yvals);
  SPHERAL_HOST_DEVICE QuadraticInterpolator() = default;

  // Initialize after construction, either with a function or tabulated values
  template<typename Func>
  void initialize(double xmin, double xmax, size_t n, const Func& f);
  void initialize(double xmin, double xmax, const std::vector<double>& yvals);

  // Comparisons
  bool operator==(const QuadraticInterpolator& rhs) const;

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
  SPHERAL_HOST_DEVICE const CoeffsType& coeffs() const;  // the fitting coefficients
  
  VVI_IMPL_DEEPCOPY(QuadraticInterpolator, mcoeffs)
  VVI_IMPL_COMPARE(QuadraticInterpolator, mN1, mXmin, mXmax, mcoeffs)
private:
  //--------------------------- Private Interface --------------------------//
  // Member data
  size_t mN1;
  double mXmin, mXmax, mXstep;
  CoeffsType mcoeffs;
};

VVI_IMPL_END

#ifdef VVI_ENABLED

class QuadraticInterpolator;

class PTR_VIEW_METACLASS_DEFAULT \
  ((QuadraticInterpolator), (QuadraticInterpolatorView), (vvimpl::QuadraticInterpolator))

#define QuadraticInterpolator__(code) PTR_VALUE_METACLASS_DECL \
  ((QuadraticInterpolator), (QuadraticInterpolatorView), (code))

class QuadraticInterpolator__(

  using CoeffsType = typename ImplType::CoeffsType;

  template<typename Func>
  QuadraticInterpolator(const double xmin,
                        const double xmax,
                        const size_t n,
                        const Func& F) 
  : VVI_VALUE_CTOR_ARGS((xmin,xmax,n,F)) {}

  VVI_VALUE_DEF_CTOR(QuadraticInterpolator)

  // Alternatively initialize from tabulated values
  void initialize(const double xmin, const double xmax,
                  const std::vector<double>& yvals)
  { VVI_IMPL_INST().initialize(xmin, xmax, yvals); }

  // Interpolate for the y value
  double operator()(const double x) const { return VVI_IMPL_INST()(x); }
  double prime(const double x) const { return VVI_IMPL_INST().prime(x); }
  double prime2(const double x) const { return VVI_IMPL_INST().prime2(x); }

  // Same as above, but use a pre-computed table position (from lowerBound)
  double operator()(const double x, const size_t i0) const { return VVI_IMPL_INST()(x, i0); }
  double prime(const double x, const size_t i0) const { return VVI_IMPL_INST().prime(x, i0); }
  double prime2(const double x, const size_t i0) const { return VVI_IMPL_INST().prime2(x, i0); }

  // Return the lower bound index in the table for the given x coordinate
  size_t lowerBound(const double x) const { return VVI_IMPL_INST().lowerBound(x); }

  // Allow read access the internal data representation
  size_t size()  const { return VVI_IMPL_INST().size(); }
  double xmin()  const { return VVI_IMPL_INST().xmin(); }
  double xmax()  const { return VVI_IMPL_INST().xmax(); }
  double xstep() const { return VVI_IMPL_INST().xstep(); }
  const vvi::vector<double>& coeffs() const { return VVI_IMPL_INST().coeffs(); }
);

#endif // VVI_ENABLED


} // namespace Spheral

#include "QuadraticInterpolatorInline.hh"

#endif
