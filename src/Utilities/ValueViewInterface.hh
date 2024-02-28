#ifndef __SPHERAL_VALUEVIEWINTERFACE__
#define __SPHERAL_VALUEVIEWINTERFACE__

#include "Field/SphArray.hh"
#include <memory>

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

template<typename T>
using SMART_PTR_TYPE = std::shared_ptr<T>;
//using SMART_PTR_TYPE = ManagedSmartPtr<T>;

template<typename view_type, typename impl_type>
class SpheralViewInterface : public SMART_PTR_TYPE<impl_type>
{
private:
  using m_ImplType = impl_type;

protected:
  using SmartPtrType = SMART_PTR_TYPE<m_ImplType>;

public:
  SPHERAL_HOST_DEVICE SmartPtrType & sptr() { return *this; }
  SPHERAL_HOST_DEVICE SmartPtrType const& sptr() const { return *this; }

  SPHERAL_HOST_DEVICE m_ImplType & sptr_data() { return *(SmartPtrType::get()); } \
  SPHERAL_HOST_DEVICE m_ImplType & sptr_data() const { return *(SmartPtrType::get()); } 

  SPHERAL_HOST_DEVICE SpheralViewInterface() = default;
  SPHERAL_HOST SpheralViewInterface(SmartPtrType&& rhs) : SmartPtrType(std::forward<SmartPtrType>(rhs)) {}
};

#define SPTR_REF ViewInterface::sptr
#define SPTR_DATA_REF ViewInterface::sptr_data

// Defines a ctor that will take a "new" Data object to create the underlying
// ManagedSmartPtr.
#define VIEW_DEFINE_ALLOC_CTOR(view_t, impl_t) \
protected: \
  view_t(impl_t* rhs) : Base(SmartPtrType(rhs, [](impl_t *p) { p->free(); } )) {}
  //view_t(impl_t* rhs) : Base(spheral::make_managedsmartptr<impl_t>(rhs)) {}

#define VIEW_DEF_CTOR(view_t) \
  view_t() = default;

#define VIEW_COPY_CTOR(view_t) \
  view_t(view_t const& rhs) = default;

#define VIEW_ASSIGNEMT_OP(view_t) \
  view_t& operator=(view_t const& rhs) = default;

#define VIEW_EQ_OP(view_t) \
  bool operator==(const view_t& rhs) const \
    { return sptr_data() == rhs.sptr_data(); }

#define VALUE_DEF_CTOR(value_t, impl_t) \
  value_t() : Base(new impl_t()) {}

#define VALUE_COPY_CTOR(value_t, impl_t) \
  value_t(value_t const& rhs) : Base(new impl_t( deepCopy( rhs.SPTR_DATA_REF() ) )) {}

#define VALUE_ASSIGNEMT_OP(value_t, impl_t) \
  value_t& operator=(value_t const& rhs) { \
    ViewType::operator=( ViewType(new impl_t( deepCopy( rhs.SPTR_DATA_REF() )))); \
    return *this; \
  }

#define VALUE_EQ_OP(value_t) \
  bool operator==(const value_t& rhs) const \
    { return compare(sptr_data(), rhs.sptr_data()); }
  

#define VALUE_TOVIEW_OP() \
  ViewType toView() { return ViewType(*this); }

#define VIEW_TYPE_ALIASES(value_t, impl) \
private: \
  using Base = SpheralViewInterface<value_t, impl>; \
  using ViewInterface = Base; \
  using ViewType = value_t; \
  using SmartPtrType = typename Base::SmartPtrType; \
  using ImplType = impl;

#define VALUE_TYPE_ALIASES(value_t, view, impl) \
private: \
  using Base = SpheralValueInterface<view, impl>; \
  using ViewInterface = SpheralViewInterface<view, impl>; \
  using ViewType = view; \
  using SmartPtrType = typename ViewInterface::SmartPtrType; \
  using ImplType = impl;

#define VALUE_TYPE_DEFAULT_ASSIGNMENT_OP(value_t) \
  value_t& operator=(value_t const& rhs) { \
    ViewType::operator=( ViewType( new ImplType(deepCopy(rhs.ViewInterface::sptr_data())) ) ); \
    return *this; \
  }

// It is assumed that the Base type is defined with from inheriting 
// SpheralViewInterface or by explicitly defining : 
//   using Base = ManagedSmartPtr<ClassType>;
#define SMART_PTR_MEMBER_ACCESSOR(type, var) \
  SPHERAL_HOST_DEVICE type & var() { return Base::get()->var; } \
  SPHERAL_HOST_DEVICE type & var() const { return Base::get()->var; }


// Interface class for Value like objects.
template<typename view_type, typename impl_type>
class SpheralValueInterface : public view_type
{
private:
  using m_ViewInterface = SpheralViewInterface<view_type, impl_type>;
  using m_ImplType = typename m_ViewInterface::SmartPtrType::element_type;
  using m_SmartPtrType = typename m_ViewInterface::SmartPtrType;

public:
  SPHERAL_HOST SpheralValueInterface(m_ImplType* rhs) : view_type((rhs)) {}
  SPHERAL_HOST SpheralValueInterface(m_SmartPtrType&& s_ptr) : view_type(std::forward<m_SmartPtrType>(s_ptr)) {}
};

//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
//  Quadratic Interpolator Example implementation
//  TODO: Remove this...
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
  VIEW_TYPE_ALIASES(QIntView, QIntData)
public:
  friend class QInt;
  using CoeffsType = typename QIntData::CoeffsType;
protected:
  VIEW_DEFINE_ALLOC_CTOR(QIntView, QIntData)
  // Interal interface for accessing the underlying members of QIntData
  SMART_PTR_MEMBER_ACCESSOR(CoeffsType, mcoeffs)

public:
  // Forward View capable methods
  SPHERAL_HOST_DEVICE double xmin() const { return sptr_data().xmin(); }
  SPHERAL_HOST_DEVICE double xmax() const { return sptr_data().xmax(); }
  SPHERAL_HOST_DEVICE CoeffsType const& coeffs() const { return sptr_data().coeffs(); }
};



class QInt : public SpheralValueInterface<QIntView, QIntData>
{
  VALUE_TYPE_ALIASES(QInt, QIntView, QIntData)
public:
  VALUE_DEF_CTOR(QInt, QIntData)
  //VALUE_COPY_CTOR(QInt, QIntData)
  //VALUE_ASSIGNEMT_OP(QInt, QIntData)
  VALUE_TOVIEW_OP()

  // Forward Value capable methods
  SPHERAL_HOST void initialize(size_t min) const { return sptr_data().initialize(min); }
  SPHERAL_HOST void editData(size_t min) const { return sptr_data().editData(min); }
  SPHERAL_HOST CoeffsType coeffs() const { return deepCopy(sptr_data().coeffs()); }
};


} // namespace Spheral

#endif // __SPHERAL_VALUEVIEWINTERFACE__
