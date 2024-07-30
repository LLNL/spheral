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
// All View classes should inherit from SpheralViewInterface. It sets up useful 
// type aliases and forwarding constructors to the underlying Data class.
//
// SpheralViewInterface children will want to use VIEW_DEFINE_ALLOC_CTOR
// in order to set up a forwarding constructor from the Value class ctor
// that passes in a "new DataObject" type.

#if defined(SPHERAL_ENABLE_SCIP)
#define SCIP_IMPL_BEGIN namespace impl {
#define SCIP_IMPL_END } // namespace impl
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

template<typename impl_type>
class SpheralViewInterface : public util::shared_ptr<impl_type>
{
private:
  using m_ImplType = impl_type;

public:
  using SmartPtrType = util::shared_ptr<m_ImplType>;
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
#if !defined(SPHERAL_ENABLE_SCIP)
  SPHERAL_HOST SpheralValueInterface(m_ImplType* rhs) : view_type((rhs)) {}
#endif
  SPHERAL_HOST SpheralValueInterface(m_SmartPtrType&& s_ptr) : view_type(std::forward<m_SmartPtrType>(s_ptr)) {}
};




#define SPTR_REF ViewInterface::sptr
#define SPTR_DATA_REF ViewInterface::sptr_data

// Defines a ctor that will take a "new" Data object to create the underlying
// ManagedSmartPtr.
#if !defined(SPHERAL_ENABLE_SCIP)
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

#define VIEW_SHALLOW_COPY(view_t) \
public: \
  void shallowCopy(view_t const& rhs) {*this = rhs;}


#define VIEW_INTERFACE_METACLASS(value_t, view_t, impl_t) \
public: \
  using ImplType = UNPACK impl_t; \
  using ViewType = UNPACK view_t; \
  using ValueType = UNPACK value_t; \
protected: \
  using Base = Spheral::SpheralViewInterface<UNPACK impl_t>; \
  using ViewInterface = Base; \
  using SmartPtrType = typename Base::SmartPtrType; \
  friend class UNPACK value_t; \
  VIEW_DEFINE_ALLOC_CTOR(view_t) \
  SPHERAL_HOST_DEVICE VIEW_DEF_CTOR(view_t) \
  VIEW_SHALLOW_COPY(UNPACK view_t) \

#define VIEW_INTERFACE_METACLASS_DECLARATION_BEGIN(value_t, view_t, impl_t) \
UNPACK view_t : public Spheral::SpheralViewInterface< UNPACK impl_t> { \
  VIEW_INTERFACE_METACLASS((UNPACK value_t), (UNPACK view_t), (UNPACK impl_t)) \
  POINTER_SYNTAX_OPERATORS()

#define VIEW_INTERFACE_METACLASS_DECLARATION_END() \
};

#define VIEW_INTERFACE_METACLASS_DECLARATION(value_t, view_t, impl_t, code) \
  VIEW_INTERFACE_METACLASS_DECLARATION_BEGIN((UNPACK value_t), (UNPACK view_t), (UNPACK impl_t)) \
  UNPACK code \
}



#define VALUE_TYPE_ALIASES(view_t) \
public: \
  using ViewType = UNPACK view_t; \
  using ValueType = typename ViewType::ValueType; \
  using ImplType = typename ViewType::ImplType; \
protected: \
  using Base = Spheral::SpheralValueInterface<ViewType>; \
  using ViewInterface = typename ViewType::ViewInterface; \
  using SmartPtrType = typename ViewType::SmartPtrType; \

#define UPCAST_CONVERSION_OP(parent_t) \
public: \
  operator parent_t() const {return parent_t(this->sptr());}

#define POINTER_SYNTAX_OPERATORS() \
public: \
  SPHERAL_HOST_DEVICE ImplType& operator*() const { return SPTR_DATA_REF(); } \
  SPHERAL_HOST_DEVICE ImplType* operator->() const { return &SPTR_DATA_REF(); }

#if !defined(SPHERAL_ENABLE_SCIP)
#define VALUE_DEF_CTOR(value_t) \
  value_t() : Base(new ImplType()) {}
#define VALUE_COPY_CTOR(value_t) \
  value_t(value_t const& rhs) : Base(new ImplType( deepCopy( rhs.SPTR_DATA_REF() ) )) {}
#define VALUE_ASSIGNEMT_OP() \
  ValueType& operator=(ValueType const& rhs) { \
    ViewType::operator=( ViewType::ViewType(new ImplType( deepCopy( rhs.SPTR_DATA_REF() )))); \
    return *this; \
  }

#else
#define VALUE_DEF_CTOR(value_t) \
  value_t() : Base(chai::make_shared<ImplType>()) {}
#define VALUE_COPY_CTOR(value_t) \
  value_t(value_t const& rhs) : Base(chai::make_shared<ImplType>( deepCopy(rhs.SPTR_DATA_REF()) )) {}
#define VALUE_ASSIGNEMT_OP() \
  ValueType& operator=(ValueType const& rhs) { \
    ViewType::operator=( ViewType::ViewType(chai::make_shared<ImplType>( deepCopy( rhs.SPTR_DATA_REF() )))); \
    return *this; \
  }

#endif


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


#define VALUE_INTERFACE_METACLASS(view_t) \
public: \
  using ViewType = UNPACK view_t; \
  using ValueType = typename ViewType::ValueType; \
  using ImplType = typename ViewType::ImplType; \
protected: \
  using Base = Spheral::SpheralValueInterface<ViewType>; \
  using ViewInterface = typename ViewType::ViewInterface; \
  using SmartPtrType = typename ViewType::SmartPtrType; \
public: \
  VALUE_TOVIEW_OP() \
  ViewType operator&() { return toView(); }\


#define VALUE_INTERFACE_METACLASS_DECLARATION_BEGIN(value_t, view_t) \
UNPACK value_t : public Spheral::SpheralValueInterface<UNPACK view_t> { \
  VALUE_INTERFACE_METACLASS((UNPACK view_t)) \


#define VALUE_INTERFACE_METACLASS_DECLARATION_END() \
};

#define VALUE_INTERFACE_METACLASS_DECLARATION(value_t, view_t, code) \
  VALUE_INTERFACE_METACLASS_DECLARATION_BEGIN((UNPACK value_t), (UNPACK view_t)) \
  UNPACK code \
}

#define DELETED_INTERFACE(type) \
public: \
  type() = delete; \
  type(type const&) = delete; \
  type& operator=(type const&) = delete; \

#define DEFAULT() ;

// It is assumed that the Base type is defined with from inheriting 
// SpheralViewInterface or by explicitly defining : 
//   using Base = ManagedSmartPtr<ClassType>;
#define SMART_PTR_MEMBER_ACCESSOR(type, var) \
  SPHERAL_HOST_DEVICE type & var() { return Base::get()->var; } \
  SPHERAL_HOST_DEVICE type & var() const { return Base::get()->var; }

} // namespace Spheral

#endif // __SPHERAL_VALUEVIEWINTERFACE__
