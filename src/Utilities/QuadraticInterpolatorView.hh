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

// This macro assumes that internal members of the XXXData type are defined as:
//   type _var; 
// It is also Assumed that the Base type is defined as: 
//   using Base = ManagedSmartPtr<ClassType>;
#define MANAGED_SMART_PTR_MEMBER_ACCESSOR(type, var) \
public: \
  SPHERAL_HOST_DEVICE type & var() { return Base::get()->_##var; } \
  SPHERAL_HOST_DEVICE type & var() const { return Base::get()->_##var; }



// Making an interface class to ensure all required methods are defined by classes.
template<typename T>
class SPHERALCopyable : public chai::CHAICopyable{
  virtual void free() = 0;
  SPHERAL_HOST_DEVICE virtual T& operator=(std::nullptr_t) = 0;
  SPHERAL_HOST_DEVICE virtual void shallowCopy(T const& rhs) = 0;
};
  


// The Data class needs to be CHAICopyable in order to trigger nested copies for Copyable 
// members within.
class QIntData : public SPHERALCopyable<QIntData>{
public:

  //SPHERAL_HOST_DEVICE QIntData() { printf("CTOR QIntData @ %p\n", this ); }
  SPHERAL_HOST_DEVICE QIntData() = default;
  SPHERAL_HOST_DEVICE QIntData(QIntData const& rhs) = default;
  //SPHERAL_HOST_DEVICE QIntData(QIntData const& rhs) : _mcoeffs(rhs._mcoeffs) {}
  SPHERAL_HOST_DEVICE QIntData& operator=(QIntData const& rhs) = default;

  //using CoeffsType = MVSmartRef<double>;
  using CoeffsType = ManagedVector<double>;

  double _mXmin, _mXmax, _mXstep;
  CoeffsType _mcoeffs;

  // We need to define a free function for usage by the ManagedSmartPtr to trigger
  // the dtor / ref_count for internally "Managed" like objects.
  void free() { _mcoeffs.free(); }
  // This type needs to be CHAICopyable so we have to define this interface for it to be usable.
  SPHERAL_HOST_DEVICE QIntData& operator=(std::nullptr_t) { _mcoeffs=nullptr; return *this; }
  SPHERAL_HOST_DEVICE void shallowCopy(QIntData const& rhs) { *this = rhs; }
};

class QIntView : public ManagedSmartPtr<QIntData> 
{
protected:
  using Base = ManagedSmartPtr<QIntData>;
  using CoeffsType = typename QIntData::CoeffsType;
  
  // This allows us to forward a new allocation of the Data we need handled into the Smart Ptr
  QIntView(QIntData* rhs) : Base(make_ManagedSmartPtr<QIntData>(rhs)) {}
  
  // Interal interface for accessing the underlying members of QIntData
  MANAGED_SMART_PTR_MEMBER_ACCESSOR(double, mXmin)
  MANAGED_SMART_PTR_MEMBER_ACCESSOR(CoeffsType, mcoeffs)

public:
  QIntView() : Base() {};
  using Base::move;

  //QIntView(QIntView const& rhs) {std::cout <<"view copyctor\n";};
  QIntView(QIntView const& rhs) = default;

  SPHERAL_HOST_DEVICE double xmin() const { return mXmin(); }
  SPHERAL_HOST_DEVICE CoeffsType const& coeffs() const { return mcoeffs(); }
};



class QInt : public QIntView
{
  using view_type = QIntView;
public:
  using view_type::CoeffsType;

  QInt() : view_type( new QIntData() ) { mcoeffs() = CoeffsType(0); }

  void initialize(size_t min)
  {
    mXmin() = min;
    mcoeffs().resize(20);
    mcoeffs().operator[](0) = 0.1;
    mcoeffs().operator[](1) = 0.2;
    mcoeffs().operator[](2) = 0.3;
    mcoeffs().operator[](19) = 0.19;
    view_type::registerTouch(chai::CPU);
  }

  view_type toView() { return view_type(*this); }
};








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
