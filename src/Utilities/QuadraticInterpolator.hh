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

class QuadraticInterpolatorImpl : SPHERALCopyable<QuadraticInterpolatorImpl> {
public:
  //--------------------------- Public Interface ---------------------------//
  //using CoeffsType = std::vector<double>;
  using CoeffsType = ManagedVector<double>;

  // Constructors, destructors
  template<typename Func>
  SPHERAL_HOST
  QuadraticInterpolatorImpl(const double xmin,
                        const double xmax,
                        const size_t n,
                        const Func& F);
  QuadraticInterpolatorImpl();
  ~QuadraticInterpolatorImpl();

  // Alternatively initialize from tabulated values
  void initialize(const double xmin, const double xmax,
                  const std::vector<double>& yvals);

  // Comparisons
  SPHERAL_HOST_DEVICE bool operator==(const QuadraticInterpolatorImpl& rhs) const;

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
  friend QuadraticInterpolatorImpl deepCopy(QuadraticInterpolatorImpl const& rhs) {
    QuadraticInterpolatorImpl result(rhs);
    result.mcoeffs = Spheral::deepCopy(rhs.mcoeffs);
    return result;
  }

  SPHERAL_HOST_DEVICE
  friend bool compare(QuadraticInterpolatorImpl const& lhs, QuadraticInterpolatorImpl const& rhs) {
    return ((lhs.mN1 == rhs.mN1) and
            (lhs.mXmin == rhs.mXmin) and
            (lhs.mXmax == rhs.mXmax) and
            (compare(lhs.mcoeffs, rhs.mcoeffs)));
  }
  void free() { mcoeffs.free(); } //Managed Vector
  SPHERAL_HOST_DEVICE QuadraticInterpolatorImpl& operator=(std::nullptr_t) { mcoeffs=nullptr; return *this; }
  SPHERAL_HOST_DEVICE void shallowCopy(QuadraticInterpolatorImpl const& rhs) { *this = rhs; }
  
  // Member data
  size_t mN1;
  double mXmin, mXmax, mXstep;
  CoeffsType mcoeffs;
};




class QuadraticInterpolatorView : 
  public SpheralViewInterface<QuadraticInterpolatorView, QuadraticInterpolatorImpl>
{
  VIEW_DEFINE_ALLOC_CTOR(QuadraticInterpolatorView, QuadraticInterpolatorImpl)

// Forward Type aliases we want to use from interfaces.
public:
  using CoeffsType = typename QuadraticInterpolatorImpl::CoeffsType;

protected:
  //SMART_PTR_MEMBER_ACCESSOR(CoeffsType, mcoeffs)

public:

  // Default Ctor
  SPHERAL_HOST_DEVICE QuadraticInterpolatorView() = default;

  // Default Copy Ctor
  SPHERAL_HOST_DEVICE QuadraticInterpolatorView(QuadraticInterpolatorView const&) = default;

  // Assignment
  SPHERAL_HOST_DEVICE QuadraticInterpolatorView& operator=(QuadraticInterpolatorView const&) = default;

  // Equivalence
  SPHERAL_HOST_DEVICE bool operator==(const QuadraticInterpolatorView& rhs) const
    { return sptr_data() == rhs.sptr_data(); }
  
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
};




class QuadraticInterpolator :
  public SpheralValueInterface<QuadraticInterpolatorView>
{
public:
  template<typename Func>
  QuadraticInterpolator(const double xmin,
                        const double xmax,
                        const size_t n,
                        const Func& F) :
    Base( new QuadraticInterpolatorImpl(xmin, xmax, n, F) ) {}

  // Default Ctor
  QuadraticInterpolator() :
    Base( new QuadraticInterpolatorImpl() ) {}

  // Copy Ctor
  QuadraticInterpolator(QuadraticInterpolator const& rhs) :
    Base( new QuadraticInterpolatorImpl(deepCopy(rhs.sptr_data())) ) {}

  // Assignment
  QuadraticInterpolator& operator=(QuadraticInterpolator const& rhs) {
    ViewType::operator=( ViewType( new QuadraticInterpolatorImpl(deepCopy(rhs.sptr_data())) ) );
    return *this;
  }
  
  // Equivalence
  bool operator==(const QuadraticInterpolator& rhs) const
    { return compare(*this, rhs); }
  

  void initialize(const double xmin, const double xmax, const std::vector<double>& yvals) 
    { sptr_data().initialize(xmin, xmax, yvals); }

  // Required interface for SpheralValueInterface classes.
  ViewType toView() { return ViewType(*this); }
};

} // namespace Spheral

#include "QuadraticInterpolatorInline.hh"

#endif
