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



//-----------------------------------------------------------------------------
// Helper classes and Macros for defining a Value/View pattern in Spheral.
//-----------------------------------------------------------------------------

// An interface class to ensure all required methods are defined by Data 
// classes that need to be automatically copied by another CHAICopyable or
// SPHERALCopyable class.
template<typename T>
class SPHERALCopyable : public chai::CHAICopyable{
  virtual void free() = 0;
  SPHERAL_HOST_DEVICE virtual T& operator=(std::nullptr_t) = 0;
  SPHERAL_HOST_DEVICE virtual void shallowCopy(T const& rhs) = 0;
};
  
// Interface class for View like objects.
//
// All View classes should inherit from SpheralViewInterface. It sets up useful 
// type aliases and forwarding constructors to the underlying Data class.
//
// SpheralViewInterface children will want to use VIEW_DEFINE_ALLOC_CTOR
// in order to set up a forwarding constructor from the Value class ctor
// that passes in a "new DataObject" type.
template<typename view_type, typename DataType>
class SpheralViewInterface : public ManagedSmartPtr<DataType>
{
protected:
  using Base = SpheralViewInterface<view_type, DataType>;
  using ViewBase = Base;

public:
  using ViewType = view_type;
  using SmartPtrType = ManagedSmartPtr<DataType>;

  SPHERAL_HOST_DEVICE DataType & sptr_data() { return *(SmartPtrType::get()); } \
  SPHERAL_HOST_DEVICE DataType & sptr_data() const { return *(SmartPtrType::get()); } 

  SPHERAL_HOST SpheralViewInterface(SmartPtrType&& rhs) : SmartPtrType(std::forward<SmartPtrType>(rhs)) {}
};

// Defines a ctor that will take a "new" Data object to create the underlying
// ManagedSmartPtr.
#define VIEW_DEFINE_ALLOC_CTOR(view_t, data_t) \
protected: \
  view_t(data_t* rhs) : ViewBase(make_ManagedSmartPtr<data_t>(rhs)) {}

// It is assumed that the Base type is defined with from inheriting 
// SpheralViewInterface or by explicitly defining : 
//   using Base = ManagedSmartPtr<ClassType>;
#define SMART_PTR_MEMBER_ACCESSOR(type, var) \
public: \
  SPHERAL_HOST_DEVICE type & var() { return ViewBase::get()->var; } \
  SPHERAL_HOST_DEVICE type & var() const { return ViewBase::get()->var; }

// Interface class for Value like objects.
template<typename view_type>
class SpheralValueInterface : public view_type
{
protected:
  using Base = SpheralValueInterface<view_type>;

private:
  using m_DataType = typename view_type::SmartPtrType::element_type;
  using m_ViewType = typename view_type::ViewType;

public:
  SPHERAL_HOST SpheralValueInterface(m_DataType* rhs) : view_type((rhs)) {}
  virtual m_ViewType toView() = 0;
};

//-----------------------------------------------------------------------------


// The Data class needs to be CHAICopyable in order to trigger nested copies for Copyable 
// members within.
class QIntData : public SPHERALCopyable<QIntData>{
public:

  SPHERAL_HOST_DEVICE QIntData() = default;
  SPHERAL_HOST_DEVICE QIntData(QIntData const& rhs) = default;
  SPHERAL_HOST_DEVICE QIntData& operator=(QIntData const& rhs) = default;

  using CoeffsType = ManagedVector<double>;

  double mXmin, mXmax, mXstep;
  CoeffsType mcoeffs;

  SPHERAL_HOST void initialize(size_t min)
  {
    mXmin = min;
    mcoeffs.resize(10);
    mcoeffs[0] = 0.1;
    mcoeffs[1] = 0.2;
    mcoeffs[2] = 0.3;
    mcoeffs[9] = 0.19;
  }

  SPHERAL_HOST void editData(size_t min)
  {
    mXmin = min;
    mcoeffs[9] = 91;
  }

  SPHERAL_HOST_DEVICE double xmin() const { return mXmin; }
  SPHERAL_HOST_DEVICE double xmax() const { return mXmax; }
  SPHERAL_HOST_DEVICE CoeffsType const& coeffs() const { return mcoeffs; }

  // Define the required interface for a SPHERALCopyable object.
  void free() { mcoeffs.free(); }
  SPHERAL_HOST_DEVICE QIntData& operator=(std::nullptr_t) { mcoeffs=nullptr; return *this; }
  SPHERAL_HOST_DEVICE void shallowCopy(QIntData const& rhs) { *this = rhs; }
};



class QIntView : public SpheralViewInterface<QIntView, QIntData>
{
  VIEW_DEFINE_ALLOC_CTOR(QIntView, QIntData)

public:
  using CoeffsType = typename QIntData::CoeffsType;
protected:
  // Interal interface for accessing the underlying members of QIntData
  SMART_PTR_MEMBER_ACCESSOR(CoeffsType, mcoeffs)

public:
  // Forward View capable methods
  SPHERAL_HOST_DEVICE double xmin() const { return sptr_data().xmin(); }
  SPHERAL_HOST_DEVICE double xmax() const { return sptr_data().xmax(); }
  SPHERAL_HOST_DEVICE CoeffsType const& coeffs() const { return sptr_data().coeffs(); }
};



class QInt : public SpheralValueInterface<QIntView>
{
public:
  QInt() : Base( new QIntData() ) {}

  // Forward Value capable methods
  SPHERAL_HOST void initialize(size_t min) const { return sptr_data().initialize(min); }
  SPHERAL_HOST void editData(size_t min) const { return sptr_data().editData(min); }
  SPHERAL_HOST_DEVICE CoeffsType coeffs() const { return deepCopy(sptr_data().coeffs()); }

  QIntView toView() { return ViewType(*this); }
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
