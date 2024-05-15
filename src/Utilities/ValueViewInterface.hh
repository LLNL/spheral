#ifndef __SPHERAL_VALUEVIEWINTERFACE__
#define __SPHERAL_VALUEVIEWINTERFACE__

#include "Utilities/SharedPtr.hh"
#include "chai/ManagedSharedPtr.hpp"
#include "Field/SphArray.hh"
#include <memory>

namespace Spheral {

// Macro tool for unpacking types with multiple template arguments.
#define UNPACK( ... ) __VA_ARGS__

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

#define USE_CMSPTR


template<typename T>
#if defined(USE_CMSPTR)
using SMART_PTR_TYPE = chai::ManagedSharedPtr<T>;
#else
using SMART_PTR_TYPE = Spheral::shared_ptr<T>;
//using SMART_PTR_TYPE = std::shared_ptr<T>;
//using SMART_PTR_TYPE = ManagedSmartPtr<T>;
#endif

template<typename view_type, typename impl_type>
class SpheralViewInterface : public SMART_PTR_TYPE<impl_type>
{
private:
  using m_ImplType = impl_type;

public:
  using SmartPtrType = SMART_PTR_TYPE<m_ImplType>;
  SPHERAL_HOST_DEVICE SmartPtrType & sptr() { return *this; }
  SPHERAL_HOST_DEVICE SmartPtrType const& sptr() const { return *this; }

  SPHERAL_HOST_DEVICE m_ImplType & sptr_data() { return *(SmartPtrType::get()); } \
  SPHERAL_HOST_DEVICE m_ImplType & sptr_data() const { return *(SmartPtrType::get()); } 

  SPHERAL_HOST_DEVICE SpheralViewInterface() = default;
  SPHERAL_HOST SpheralViewInterface(SmartPtrType&& rhs) : SmartPtrType(std::forward<SmartPtrType>(rhs)) {}
};


// Interface class for Value like objects.
template<typename view_type>
class SpheralValueInterface : public view_type
{
private:
  using m_ImplType = typename view_type::ImplType;
  using m_SmartPtrType = typename view_type::SmartPtrType;

protected:
  using ViewType = typename view_type::ViewType;

public:
#if !defined(USE_CMSPTR)
  SPHERAL_HOST SpheralValueInterface(m_ImplType* rhs) : view_type((rhs)) {}
#endif
  SPHERAL_HOST SpheralValueInterface(m_SmartPtrType&& s_ptr) : view_type(std::forward<m_SmartPtrType>(s_ptr)) {}
};




#define SPTR_REF ViewInterface::sptr
#define SPTR_DATA_REF ViewInterface::sptr_data

#define VIEW_TYPE_ALIASES(value_t, view_t, impl_t) \
public: \
  using ImplType = UNPACK impl_t; \
  using ViewType = UNPACK view_t; \
  using ValueType = UNPACK value_t; \
protected: \
  using Base = Spheral::SpheralViewInterface<UNPACK view_t, UNPACK impl_t>; \
  using ViewInterface = Base; \
  using SmartPtrType = typename Base::SmartPtrType; \
  friend class UNPACK value_t;

#define VALUE_TYPE_ALIASES(view_t) \
public: \
  using ViewType = UNPACK view_t; \
  using ValueType = typename ViewType::ValueType; \
  using ImplType = typename ViewType::ImplType; \
protected: \
  using Base = Spheral::SpheralValueInterface<ViewType>; \
  using ViewInterface = typename ViewType::ViewInterface; \
  using SmartPtrType = typename ViewType::SmartPtrType; \

// Defines a ctor that will take a "new" Data object to create the underlying
// ManagedSmartPtr.
#if !defined(USE_CMSPTR)
#define VIEW_DEFINE_ALLOC_CTOR(view_t) \
protected: \
  view_t(ImplType* rhs) : Base(SmartPtrType(rhs, [](ImplType *p) { p->free(); } )) {}
  //view_t(ImplType* rhs) : Base(SmartPtrType(rhs, [](ImplType *p) { p->free(); } )) {}
#else
#define VIEW_DEFINE_ALLOC_CTOR(view_t) \
public: \
  view_t(SmartPtrType&& rhs) : Base(std::forward<SmartPtrType>(rhs)) {}
#endif

#define VIEW_DEF_CTOR(view_t) \
  view_t() = default;

#define VIEW_COPY_CTOR(view_t) \
  view_t(view_t const& rhs) = default;

#define VIEW_ASSIGNEMT_OP() \
  SPHERAL_HOST_DEVICE ViewType& operator=(ViewType const& rhs) = default;

#define VIEW_EQ_OP() \
  bool operator==(const ViewType& rhs) const \
    { return sptr_data() == rhs.sptr_data(); }

#if !defined(USE_CMSPTR)
#define VALUE_DEF_CTOR(value_t) \
  value_t() : Base(new ImplType()) {}
#else
#define VALUE_DEF_CTOR(value_t) \
  value_t() : Base(chai::make_shared<ImplType>()) {}
#endif

#define VALUE_COPY_CTOR(value_t) \
  value_t(value_t const& rhs) : Base(new ImplType( deepCopy( rhs.SPTR_DATA_REF() ) )) {}

#define VALUE_ASSIGNEMT_OP() \
  ValueType& operator=(ValueType const& rhs) { \
    ViewType::operator=( ViewType::ViewType(new ImplType( deepCopy( rhs.SPTR_DATA_REF() )))); \
    return *this; \
  }

#define VALUE_EQ_OP() \
  bool operator==(const ValueType& rhs) const \
    { return compare(sptr_data(), rhs.sptr_data()); }

#define VALUE_TOVIEW_OP() \
  ViewType toView() { return ViewType(*this); }

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


class QInt;

class QIntView : public SpheralViewInterface<QIntView, QIntData>
{
  VIEW_TYPE_ALIASES((QInt), (QIntView), (QIntData))
public:
  friend class QInt;
  using CoeffsType = typename QIntData::CoeffsType;
protected:
  VIEW_DEFINE_ALLOC_CTOR(QIntView)
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
  VALUE_TYPE_ALIASES((QIntView))
public:
  VALUE_DEF_CTOR(QInt)
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
