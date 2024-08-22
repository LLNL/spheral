#ifndef __SPHERAL_VALUEVIEWINTERFACE__
#define __SPHERAL_VALUEVIEWINTERFACE__

//#include "Utilities/SharedPtr.hh"
#include "ManagedVector.hh"
#include "config.hh"
#include "chai/ManagedSharedPtr.hpp"
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
class SPHERALCopyable : public chai::CHAICopyable{};
  
// Interface class for View like objects.
//
// All View classes should inherit from ViewInterface. It sets up useful 
// type aliases and forwarding constructors to the underlying Data class.
//
// ViewInterface children will want to use VIEW_DEFINE_ALLOC_CTOR
// in order to set up a forwarding constructor from the Value class ctor
// that passes in a "new DataObject" type.

#if defined(SPHERAL_ENABLE_SCIP)
#define SCIP_IMPL_BEGIN namespace vvi { namespace impl {
#define SCIP_IMPL_END } } // namespace vvi::impl
#else
#define SCIP_IMPL_BEGIN
#define SCIP_IMPL_END
#endif

namespace util {

template<typename T>
#if defined(SPHERAL_ENABLE_SCIP)
using shared_ptr = chai::ManagedSharedPtr<T>;
#else
using shared_ptr = std::shared_ptr<T>;
#endif

template<typename T>
#if defined(SPHERAL_ENABLE_SCIP)
using vector = ManagedVector<T>;
#else
using vector = std::vector<T>;
#endif

}// namespace util

} // namespace Spheral


namespace vvi {

template<typename impl_type>
class ViewInterface : public ::Spheral::util::shared_ptr<impl_type>
{
private:
  using m_ImplType = impl_type;

public:
  using SmartPtrType = ::Spheral::util::shared_ptr<m_ImplType>;
  SPHERAL_HOST_DEVICE SmartPtrType & sptr() { return *this; }
  SPHERAL_HOST_DEVICE SmartPtrType const& sptr() const { return *this; }

  SPHERAL_HOST_DEVICE m_ImplType & sptr_data() { return *(SmartPtrType::get()); } \
  SPHERAL_HOST_DEVICE m_ImplType & sptr_data() const { return *(SmartPtrType::get()); } 

  SPHERAL_HOST_DEVICE ViewInterface() = default;
  SPHERAL_HOST ViewInterface(SmartPtrType&& rhs) : SmartPtrType(std::forward<SmartPtrType>(rhs)) {}
};


// Interface class for Value like objects.
template<typename view_type>
class ValueInterface : public view_type
{
private:
  using m_ImplType = typename view_type::ImplType;
  using m_SmartPtrType = typename view_type::SmartPtrType;

protected:
  using ViewType = typename view_type::ViewType;

public:
#if !defined(SPHERAL_ENABLE_SCIP)
  SPHERAL_HOST ValueInterface(m_ImplType* rhs) : view_type((rhs)) {}
#endif
  SPHERAL_HOST ValueInterface(m_SmartPtrType&& s_ptr) : view_type(std::forward<m_SmartPtrType>(s_ptr)) {}
};

} // namespace vvi

#define vvi_impl_SPTR_DATA_REF ViewInterfaceType::sptr_data






// Defines a ctor that will take a "new" Data object to create the underlying
// ManagedSmartPtr.
#if !defined(SPHERAL_ENABLE_SCIP)
#define vvi_impl_VIEW_DEFINE_ALLOC_CTOR(view_t) \
protected: \
  view_t(ImplType* rhs) : Base(SmartPtrType(rhs, [](ImplType *p) { p->free(); } )) {}
  //view_t(ImplType* rhs) : Base(SmartPtrType(rhs, [](ImplType *p) { p->free(); } )) {}
#else
#define vvi_impl_VIEW_DEFINE_ALLOC_CTOR(view_t) \
public: \
  view_t(SmartPtrType&& rhs) : Base(std::forward<SmartPtrType>(rhs)) {}
#endif

#define vvi_impl_VIEW_SHALLOW_COPY(view_t) \
public: \
  void shallowCopy(view_t const& rhs) {*this = rhs;}





// ----------------------------------------------------------------------------
// Behavioral macros : These are used to allow for pointer like behavior from
// interface classes
// ----------------------------------------------------------------------------

#define POINTER_SYNTAX_OPERATORS() \
public: \
  SPHERAL_HOST_DEVICE ImplType& operator*() const { return vvi_impl_SPTR_DATA_REF(); } \
  SPHERAL_HOST_DEVICE ImplType* operator->() const { return &vvi_impl_SPTR_DATA_REF(); }


#define REF_OPERATOR() \
public: \
  ViewType operator&() { return toView(); }


#define VVI_UPCAST_CONVERSION_OP(parent_t) \
public: \
  operator parent_t() const {return parent_t(this->sptr());}


#define VVI_DELETED_INTERFACE(type) \
public: \
  type() = delete; \
  type(type const&) = delete; \
  type& operator=(type const&) = delete;

#define VVI_DEFAULT() ;

#define VVI_MEMBER_ACCESSOR(type, var) \
  SPHERAL_HOST_DEVICE type & var() { return Base::get()->var; } \
  SPHERAL_HOST_DEVICE type & var() const { return Base::get()->var; }

#define VVI_IMPL_INST this->sptr_data()

// ----------------------------------------------------------------------------
// VIEW class definition macros
// ----------------------------------------------------------------------------

#define VVI_VIEW_DEF_CTOR(view_t) \
  view_t() = default;

#define VVI_VIEW_COPY_CTOR(view_t) \
  view_t(view_t const& rhs) = default;

#define VVI_VIEW_ASSIGNEMT_OP() \
  SPHERAL_HOST_DEVICE ViewType& operator=(ViewType const& rhs) = default;

#define VVI_VIEW_EQ_OP() \
  bool operator==(const ViewType& rhs) const \
    { return sptr_data() == rhs.sptr_data(); }




// ----------------------------------------------------------------------------
// VIEW class declaration macros
// ----------------------------------------------------------------------------

#define VIEW_INTERFACE_METACLASS(value_t, view_t, impl_t) \
public: \
  using ImplType = UNPACK impl_t; \
  using ViewType = UNPACK view_t; \
  using ValueType = UNPACK value_t; \
protected: \
  using Base = ::vvi::ViewInterface<UNPACK impl_t>; \
  using ViewInterfaceType = Base; \
  using SmartPtrType = typename Base::SmartPtrType; \
  friend class UNPACK value_t; \
  vvi_impl_VIEW_DEFINE_ALLOC_CTOR(view_t) \
private:

#define VIEW_METACLASS_DECL_BEGIN(value_t, view_t, impl_t) \
UNPACK view_t : public ::vvi::ViewInterface< UNPACK impl_t> { \
  VIEW_INTERFACE_METACLASS((UNPACK value_t), (UNPACK view_t), (UNPACK impl_t))

#define VIEW_METACLASS_DECL_END() \
};

#define PTR_VIEW_METACLASS_DECL(value_t, view_t, impl_t, code) \
  VIEW_METACLASS_DECL_BEGIN((UNPACK value_t), (UNPACK view_t), (UNPACK impl_t)) \
  POINTER_SYNTAX_OPERATORS() \
  UNPACK code \
  VIEW_METACLASS_DECL_END()

#define REF_VIEW_METACLASS_DECL(value_t, view_t, impl_t, code) \
  VIEW_METACLASS_DECL_BEGIN((UNPACK value_t), (UNPACK view_t), (UNPACK impl_t)) \
  UNPACK code \
  VIEW_METACLASS_DECL_END()








// ----------------------------------------------------------------------------
// VALUE class definition macros
// ----------------------------------------------------------------------------

#define VVI_CTOR(value_t, params, args) \
  value_t(UNPACK params) : Base(chai::make_shared<ImplType>(UNPACK args)) {}

#define VVI_CTOR_ARGS(args) \
  Base(chai::make_shared<ImplType>(UNPACK args))

#define VVI_DTOR(value_t, code) \
  ~value_t() { UNPACK code }

#define VVI_VALUE_DEF_CTOR(value_t) \
  value_t() : Base(chai::make_shared<ImplType>()) {}

#define VVI_VALUE_COPY_CTOR(value_t) \
  value_t(value_t const& rhs) : Base(chai::make_shared<ImplType>( deepCopy(rhs.vvi_impl_SPTR_DATA_REF()) )) {}

#define VVI_VALUE_ASSIGNEMT_OP() \
  ValueType& operator=(ValueType const& rhs) { \
    ViewType::operator=( ViewType::ViewType(chai::make_shared<ImplType>( deepCopy( rhs.vvi_impl_SPTR_DATA_REF() )))); \
    return *this; \
  }

#define VVI_VALUE_EQ_OP() \
  bool operator==(const ValueType& rhs) const \
    { return compare(sptr_data(), rhs.sptr_data()); }

#define VVI_VALUE_TOVIEW_OP() \
  ViewType toView() { return ViewType(*this); }

#define VVI_VALUE_TYPE_DEFAULT_ASSIGNMENT_OP(value_t) \
  value_t& operator=(value_t const& rhs) { \
    ViewType::operator=( ViewType( new ImplType(deepCopy(rhs.ViewInterfaceType::sptr_data())) ) ); \
    return *this; \
  }


// ----------------------------------------------------------------------------
// VALUE class declaration macros
// ----------------------------------------------------------------------------

#define VALUE_INTERFACE_METACLASS(view_t) \
public: \
  using ViewType = UNPACK view_t; \
  using ValueType = typename ViewType::ValueType; \
  using ImplType = typename ViewType::ImplType; \
protected: \
  using Base = ::vvi::ValueInterface<ViewType>; \
  using ViewInterfaceType = typename ViewType::ViewInterfaceType; \
  using SmartPtrType = typename ViewType::SmartPtrType; \
public: \
  VVI_VALUE_TOVIEW_OP() \
private:


#define VALUE_METACLASS_DECL_BEGIN(value_t, view_t) \
UNPACK value_t : public ::vvi::ValueInterface<UNPACK view_t> { \
  VALUE_INTERFACE_METACLASS((UNPACK view_t)) \


#define VALUE_METACLASS_DECL_END() \
};

#define PTR_VALUE_METACLASS_DECL(value_t, view_t, code) \
  VALUE_METACLASS_DECL_BEGIN((UNPACK value_t), (UNPACK view_t)) \
  REF_OPERATOR() \
  UNPACK code \
  VALUE_METACLASS_DECL_END()


#define REF_VALUE_METACLASS_DECL(value_t, view_t, code) \
  VALUE_METACLASS_DECL_BEGIN((UNPACK value_t), (UNPACK view_t)) \
  UNPACK code \
  VALUE_METACLASS_DECL_END()





#endif // __SPHERAL_VALUEVIEWINTERFACE__
