#ifndef __SPHERAL_VALUEVIEWINTERFACE__
#define __SPHERAL_VALUEVIEWINTERFACE__

#include "config.hh"
#include "ManagedVector.hh"
#include "chai/ManagedSharedPtr.hpp"
#include "ValueViewInterfaceImpl.hh"

#include <memory>


namespace Spheral {

// An interface cass to ensure all required methods are defined by Data 
// classes that need to be automatically copied by another CHAICopyable or
// SPHERALCopyable class.
using SPHERALCopyable = chai::CHAICopyable;
  
} // namespace Spheral

// ----------------------------------------------------------------------------
// IMPL class declaration macros
// ----------------------------------------------------------------------------

#if defined(VVI_ENABLED)
#define VVI_IMPL_BEGIN namespace vvimpl {
#define VVI_IMPL_END } // namespace vvimpl
#else
#define VVI_IMPL_BEGIN
#define VVI_IMPL_END
#endif // defined(VVI_ENABLED)

#if defined(VVI_ENABLED)

#define VVI_IMPL_DEEPCOPY(...) VVI_IMPL_DEEPCOPY__( __VA_ARGS__ , void)
#define VVI_IMPL_COMPARE(...) VVI_IMPL_COMPARE__( __VA_ARGS__, void) 

#else

#define VVI_IMPL_DEEPCOPY(...)
#define VVI_IMPL_COMPARE(...)

#endif

// ----------------------------------------------------------------------------
// VALUE class declaration macros
// ----------------------------------------------------------------------------

#define PTR_VALUE_METACLASS_DECL(value_t, view_t, code) \
  VALUE_METACLASS_DECL_BEGIN((UNPACK value_t), (UNPACK view_t)) \
  REF_OPERATOR() \
  UNPACK code \
  VALUE_METACLASS_DECL_END()


#define REF_VALUE_METACLASS_DECL(value_t, view_t, code) \
  VALUE_METACLASS_DECL_BEGIN((UNPACK value_t), (UNPACK view_t)) \
  UNPACK code \
  VALUE_METACLASS_DECL_END()

// ----------------------------------------------------------------------------
// VIEW class declaration macros
// ----------------------------------------------------------------------------

#define PTR_VIEW_METACLASS_DECL(value_t, view_t, impl_t, code) \
  VIEW_METACLASS_DECL_BEGIN((UNPACK value_t), (UNPACK view_t), (UNPACK impl_t)) \
  POINTER_SYNTAX_OPERATORS() \
  UNPACK code \
  VIEW_METACLASS_DECL_END()

#define REF_VIEW_METACLASS_DECL(value_t, view_t, impl_t, code) \
  VIEW_METACLASS_DECL_BEGIN((UNPACK value_t), (UNPACK view_t), (UNPACK impl_t)) \
  UNPACK code \
  VIEW_METACLASS_DECL_END()

#define PTR_VIEW_METACLASS_DEFAULT(value_t, view_t, impl_t) \
  PTR_VIEW_METACLASS_DECL( value_t, view_t, impl_t, (  ) );

#define PTR_VALUE_METACLASS_DEFAULT(value_t, view_t, impl_t) \
  PTR_VALUE_METACLASS_DECL( value_t, view_t, impl_t, (  ) );

#define PTR_VALUE_METACLASS_DELETED(value_t, view_t, impl_t) \
  PTR_VALUE_METACLASS_DECL( value_t, view_t, ( VVI_DELETED_INTERFACE(UNPACK value_t) ) );

// ----------------------------------------------------------------------------
// VALUE class definition macros
// ----------------------------------------------------------------------------

#define VVI_VALUE_CTOR_ARGS(args) \
  Base(VVI_MAKE_SHARED<ImplType>(UNPACK args))

#define VVI_VALUE_DEF_CTOR(value_t) \
  value_t() : Base() {}

#define VVI_VALUE_COPY_CTOR(value_t) \
  value_t(value_t const& rhs) : Base(rhs) {}

/*
#define VVI_VALUE_ASSIGNEMT_OP() \
  ValueType& operator=(ValueType const& rhs) { \
    Base::operator=(rhs); \
    return *this; \
  }

#define VVI_VALUE_EQ_OP() \
  bool operator==(const ValueType& rhs) const \
    { return compare(VVI_SPTR_DATA_REF__(), rhs.VVI_SPTR_DATA_REF__()); }

#define VVI_VALUE_TYPE_DEFAULT_ASSIGNMENT_OP(value_t) \
  value_t& operator=(value_t const& rhs) { \
    ViewType::operator=( ViewType( new ImplType(deepCopy(rhs.ViewInterfaceType::sptr_data())) ) ); \
    return *this; \
  }
*/

// ----------------------------------------------------------------------------
// VIEW class definition macros
// ----------------------------------------------------------------------------

/*
#define VVI_VIEW_COPY_CTOR(view_t) \
  view_t(view_t const& rhs) = default;

#define VVI_VIEW_ASSIGNEMT_OP() \
  SPHERAL_HOST_DEVICE ViewType& operator=(ViewType const& rhs) = default;

#define VVI_VIEW_EQ_OP() \
  bool operator==(const ViewType& rhs) const \
    { return sptr_data() == rhs.sptr_data(); }
*/

// ----------------------------------------------------------------------------
// Special case macros
// ----------------------------------------------------------------------------

// This is used for allowing implicit upcast conversions to a Base Class
#define VVI_UPCAST_CONVERSION_OP(parent_t) \
public: \
  operator parent_t() const {return parent_t(this->sptr());}

// This is used when declaring a Value interface of an abstract class.
#define VVI_DELETED_INTERFACE(type) \
public: \
  type() = delete; \
  type(type const&) = delete; \
  type& operator=(type const&) = delete;

/*
// Used for allowing default class behavior from the compiler. Typically this
// will be used for invoking a basic view interface when using poiter syntax.
#define VVI_VIEW_DEFAULT(view_t) \
  VVI_VIEW_DEF_CTOR(view_t) \
  VVI_VIEW_COPY_CTOR(view_t) \
  VVI_VIEW_ASSIGNEMT_OP()
*/

// Creates a function that allows interface objects to access member variables.
// through a function of the same name as the variable itself.
#define VVI_MEMBER_ACCESSOR(type, var) \
  SPHERAL_HOST_DEVICE type & var() { return Base::get()->var; } \
  SPHERAL_HOST_DEVICE type & var() const { return Base::get()->var; }


// Grab the actual insatance of the implementation class the interface is accessing.
// Used when calling out to functions in the interface definitions.
//#define VVI_IMPL_INST this->sptr_data()
#define VVI_IMPL_INST VVI_SPTR_DATA_REF__



#endif // __SPHERAL_VALUEVIEWINTERFACE__
