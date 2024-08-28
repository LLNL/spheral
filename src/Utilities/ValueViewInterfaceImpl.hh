#ifndef __SPHERAL_VALUEVIEWINTERFACEIMPL__
#define __SPHERAL_VALUEVIEWINTERFACEIMPL__

// Macro tool for unpacking types with multiple template arguments.
#define UNPACK( ... ) __VA_ARGS__



// ----------------------------------------------------------------------------
// Class definitions for Value and View Intreface class structure.
// ----------------------------------------------------------------------------

namespace vvi {

template<typename impl_type>
class ViewInterface : public vvi::shared_ptr<impl_type>
{
private:
  using m_ImplType = impl_type;

public:
  using SmartPtrType = vvi::shared_ptr<m_ImplType>;
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
#if !defined(SPHERAL_ENABLE_VVI)
  SPHERAL_HOST ValueInterface(m_ImplType* rhs) : view_type((rhs)) {}
#endif
  SPHERAL_HOST ValueInterface(m_SmartPtrType&& s_ptr) : view_type(std::forward<m_SmartPtrType>(s_ptr)) {}
};

} // namespace vvi


// Internale macro for accessing data
#define VVI_SPTR_DATA_REF__ ViewInterfaceType::sptr_data


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



// ----------------------------------------------------------------------------
// VIEW class declaration macros
// ----------------------------------------------------------------------------

// Defines a ctor that will take a "new" Data object to create the underlying
// ManagedSmartPtr.
#if !defined(SPHERAL_ENABLE_VVI)
#define VVI_VIEW_DEFINE_ALLOC_CTOR__(view_t) \
protected: \
  view_t(ImplType* rhs) : Base(SmartPtrType(rhs, [](ImplType *p) { p->free(); } )) {}
  //view_t(ImplType* rhs) : Base(SmartPtrType(rhs, [](ImplType *p) { p->free(); } )) {}
#else
#define VVI_VIEW_DEFINE_ALLOC_CTOR__(view_t) \
public: \
  view_t(SmartPtrType&& rhs) : Base(std::forward<SmartPtrType>(rhs)) {}
#endif

// TODO: Is this still relevant?
//#define VVI_VIEW_SHALLOW_COPY__(view_t) \
//public: \
//  void shallowCopy(view_t const& rhs) {*this = rhs;}

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
  VVI_VIEW_DEFINE_ALLOC_CTOR__(view_t) \
private:

#define VIEW_METACLASS_DECL_BEGIN(value_t, view_t, impl_t) \
UNPACK view_t : public ::vvi::ViewInterface< UNPACK impl_t> { \
  VIEW_INTERFACE_METACLASS((UNPACK value_t), (UNPACK view_t), (UNPACK impl_t))

#define VIEW_METACLASS_DECL_END() \
};


// ----------------------------------------------------------------------------
// Behavioral macros : These are used to allow for pointer like behavior from
// interface classes
// ----------------------------------------------------------------------------

#define POINTER_SYNTAX_OPERATORS() \
public: \
  SPHERAL_HOST_DEVICE ImplType& operator*() const { return VVI_SPTR_DATA_REF__(); } \
  SPHERAL_HOST_DEVICE ImplType* operator->() const { return &VVI_SPTR_DATA_REF__(); }


#define REF_OPERATOR() \
public: \
  ViewType operator&() { return toView(); }


#endif // __SPHERAL_VALUEVIEWINTERFACEIMPL__

