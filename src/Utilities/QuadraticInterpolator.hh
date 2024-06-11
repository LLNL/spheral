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
#include "Utilities/ValueViewInterface.hh"
#include "ValueViewInterface.hh"

namespace Spheral {

namespace impl {

class QuadraticInterpolator : SPHERALCopyable<QuadraticInterpolator> {
public:
  //--------------------------- Public Interface ---------------------------//
  //using CoeffsType = std::vector<double>;
  using CoeffsType = Spheral::ManagedVector<double>;

  // Constructors, destructors
  template<typename Func>
  SPHERAL_HOST
  QuadraticInterpolator(const double xmin,
                        const double xmax,
                        const size_t n,
                        const Func& F);
  SPHERAL_HOST_DEVICE
  QuadraticInterpolator() = default;
  SPHERAL_HOST_DEVICE
  ~QuadraticInterpolator() = default;

  // Alternatively initialize from tabulated values
  void initialize(const double xmin, const double xmax,
                  const std::vector<double>& yvals);

  // Comparisons
  SPHERAL_HOST_DEVICE bool operator==(const QuadraticInterpolator& rhs) const;

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
                                              
  // Define the required interface for a SPHERALCopyable object.
  friend QuadraticInterpolator deepCopy(QuadraticInterpolator const& rhs) {
    QuadraticInterpolator result(rhs);
    result.mcoeffs = Spheral::deepCopy(rhs.mcoeffs);
    return result;
  }

  SPHERAL_HOST_DEVICE
  friend bool compare(QuadraticInterpolator const& lhs, QuadraticInterpolator const& rhs) {
    return ((lhs.mN1 == rhs.mN1) and
            (lhs.mXmin == rhs.mXmin) and
            (lhs.mXmax == rhs.mXmax) and
            (compare(lhs.mcoeffs, rhs.mcoeffs)));
  }
  void free() { mcoeffs.free(); } //Managed Vector
  SPHERAL_HOST_DEVICE QuadraticInterpolator& operator=(std::nullptr_t) { mcoeffs=nullptr; return *this; }
  SPHERAL_HOST_DEVICE void shallowCopy(QuadraticInterpolator const& rhs) { *this = rhs; }
  
  // Member data
  size_t mN1;
  double mXmin, mXmax, mXstep;
  CoeffsType mcoeffs;
};

} // namespace impl


class QuadraticInterpolator;

#define QuadraticInterpolatorView__(code) \
  VIEW_INTERFACE_METACLASS_DECLARATION((QuadraticInterpolator), \
                                       (QuadraticInterpolatorView), \
                                       (impl::QuadraticInterpolator), \
                                       (code))
#define QuadraticInterpolator__(code) \
  VALUE_INTERFACE_METACLASS_DECLARATION((QuadraticInterpolator), \
                                        (QuadraticInterpolatorView), \
                                        (code))


class QuadraticInterpolatorView__(
public:
  SPHERAL_HOST_DEVICE VIEW_COPY_CTOR(QuadraticInterpolatorView)
  SPHERAL_HOST_DEVICE VIEW_ASSIGNEMT_OP()
  SPHERAL_HOST_DEVICE VIEW_EQ_OP()
);


class QuadraticInterpolator__(
public:
  VALUE_DEF_CTOR(QuadraticInterpolator)
  VALUE_COPY_CTOR(QuadraticInterpolator)
  VALUE_ASSIGNEMT_OP()
  VALUE_EQ_OP()

  using CoeffsType = typename ImplType::CoeffsType;

  template<typename Func>
  QuadraticInterpolator(const double xmin,
                        const double xmax,
                        const size_t n,
                        const Func& F) :
    Base( chai::make_shared<ImplType>(xmin, xmax, n, F) ) {}

  void initialize(const double xmin, const double xmax, const std::vector<double>& yvals) 
    { sptr_data().initialize(xmin, xmax, yvals); }

  SPHERAL_HOST_DEVICE double operator()(const double x) const { return sptr_data()(x); }
  SPHERAL_HOST_DEVICE double prime(const double x) const { return sptr_data().prime(x); }
  SPHERAL_HOST_DEVICE double prime2(const double x) const { return sptr_data().prime2(x); }

  SPHERAL_HOST_DEVICE double operator()(const double x, const size_t i0) const { return sptr_data()(x, i0); }
  SPHERAL_HOST_DEVICE double prime(const double x, const size_t i0) const { return sptr_data().prime(x, i0); }
  SPHERAL_HOST_DEVICE double prime2(const double x, const size_t i0) const{ return sptr_data().prime2(x, i0); }

  SPHERAL_HOST_DEVICE size_t lowerBound(const double x) const { return sptr_data().lowerBound(x); }

  SPHERAL_HOST_DEVICE double size() const { return sptr_data().size(); }
  SPHERAL_HOST_DEVICE double xmin() const { return sptr_data().xmin(); }
  SPHERAL_HOST_DEVICE double xmax() const { return sptr_data().xmax(); }
  SPHERAL_HOST_DEVICE double xstep() const { return sptr_data().xstep(); }
  SPHERAL_HOST_DEVICE CoeffsType const& coeffs() const { return sptr_data().coeffs(); }
);

} // namespace Spheral

#include "QuadraticInterpolatorInline.hh"

#endif
