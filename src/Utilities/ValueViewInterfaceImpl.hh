#ifndef __SPHERAL_VALUEVIEWINTERFACEIMPL__
#define __SPHERAL_VALUEVIEWINTERFACEIMPL__

// Macro tool for unpacking types with multiple template arguments.
#define UNPACK( ... ) __VA_ARGS__

// ----------------------------------------------------------------------------
// Class definitions for Value and View Intreface class structure.
// ----------------------------------------------------------------------------

namespace vvi {

template<typename T>
#if defined(VVI_ENABLED)
using shared_ptr = chai::ManagedSharedPtr<T>;
#define VVI_MAKE_SHARED chai::make_shared
#else
using shared_ptr = std::shared_ptr<T>;
#define VVI_MAKE_SHARED std::make_shared
#endif

template<typename T>
#if defined(VVI_ENABLED)
using vector = ::Spheral::ManagedVector<T>;
//using vector = std::vector<T>;
#else
using vector = std::vector<T>;
#endif

namespace detail {

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

    ValueInterface() : ValueInterface(VVI_MAKE_SHARED<m_ImplType>()) {}
    //ValueInterface() : ValueInterface(chai::make_shared<m_ImplType>()) {}

  public:
  #if !defined(SPHERAL_ENABLE_VVI)
    ValueInterface(m_ImplType* rhs) : view_type((rhs)) {}
  #endif
    ValueInterface(m_SmartPtrType&& s_ptr) : view_type(std::forward<m_SmartPtrType>(s_ptr)) {std::cout << "fwd Ctor\n";}
    ValueInterface(ValueInterface const& rhs) : ValueInterface(VVI_MAKE_SHARED<m_ImplType>( deepCopy(rhs.sptr_data()) )) {std::cout << "Copy Ctor\n";}
    ValueInterface& operator=(ValueInterface const& rhs) { ViewType::operator=( VVI_MAKE_SHARED<m_ImplType>( deepCopy( rhs.sptr_data() ) ) ); std::cout << "Ass Op/n"; return *this; }
    bool operator==(ValueInterface const& rhs) const { return compare(this->sptr_data(), rhs.sptr_data()); }
  };


  template<typename T>
  bool compare(T const& lhs, T const& rhs) 
  { return lhs == rhs; }

  template<typename elem>
  bool compare(Spheral::ManagedVector<elem> const& lhs, Spheral::ManagedVector<elem> const& rhs)
  { return ::Spheral::compare(lhs, rhs); }

}

} // namespace vvi


// Internale macro for accessing data
#define VVI_SPTR_DATA_REF__ ViewInterfaceType::sptr_data

// ----------------------------------------------------------------------------
// IMPL class declaration macros
// ----------------------------------------------------------------------------

// Macro tool for creating macro functions with an unknown number of arguments (Max 9)
#define MKFNS(fn,...) MKFN_N(fn,##__VA_ARGS__, 9,8,7,6,5,4,3,2,1,0, void)(__VA_ARGS__)
#define MKFN_N(fn, n0, n1, n2, n3, n4, n5, n6, n7, n8, n9, n, ...) fn##n


// DeepCopy macro helpers

#define VVI_IMPL_DEEPCOPY__(impl_t, ...) \
  friend impl_t deepCopy(impl_t const& rhs) { \
    impl_t result(rhs); \
    VVI_IDC_EXPAND__( __VA_ARGS__ ) \
    return result; \
  }

#define VIDCH__(arg) result.arg = deepCopy(rhs.arg);

#define VVI_IDC_EXPAND__(...) MKFNS(VVI_IDC_EXPAND__,##__VA_ARGS__)
#define VVI_IDC_EXPAND__0(void) 
#define VVI_IDC_EXPAND__1(arg1, void) \
  VIDCH__(arg1)
#define VVI_IDC_EXPAND__2(arg1, arg2, void) \
  VIDCH__(arg1) VIDCH__(arg2)
#define VVI_IDC_EXPAND__3(arg1, arg2, arg3, void) \
  VIDCH__(arg1) VIDCH__(arg2) VIDCH__(arg3)
#define VVI_IDC_EXPAND__4(arg1, arg2, arg3, arg4, void) \
  VIDCH__(arg1) VIDCH__(arg2) VIDCH__(arg3) \
  VIDCH__(arg4) 
#define VVI_IDC_EXPAND__5(arg1, arg2, arg3, arg4, arg5, void) \
  VIDCH__(arg1) VIDCH__(arg2) VIDCH__(arg3) \
  VIDCH__(arg4) VIDCH__(arg5)
#define VVI_IDC_EXPAND__6(arg1, arg2, arg3, arg4, arg5, arg6, void) \
  VIDCH__(arg1) VIDCH__(arg2) VIDCH__(arg3) \
  VIDCH__(arg4) VIDCH__(arg5) VIDCH__(arg6)
#define VVI_IDC_EXPAND__7(arg1, arg2, arg3, arg4, arg5, arg6, arg7, void) \
  VIDCH__(arg1) VIDCH__(arg2) VIDCH__(arg3) \
  VIDCH__(arg4) VIDCH__(arg5) VIDCH__(arg6) \
  VIDCH__(arg7) 
#define VVI_IDC_EXPAND__8(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, void) \
  VIDCH__(arg1) VIDCH__(arg2) VIDCH__(arg3) \
  VIDCH__(arg4) VIDCH__(arg5) VIDCH__(arg6) \
  VIDCH__(arg7) VIDCH__(arg8)
#define VVI_IDC_EXPAND__9(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, void) \
  VIDCH__(arg1) VIDCH__(arg2) VIDCH__(arg3) \
  VIDCH__(arg4) VIDCH__(arg5) VIDCH__(arg6) \
  VIDCH__(arg7) VIDCH__(arg8) VIDCH__(arg9)

// Compare macro helpers

#define VVI_IMPL_COMPARE__(impl_t, ...) \
  friend bool compare(impl_t const& lhs, impl_t const& rhs) { \
    return VVI_ICOM_EXPAND__( __VA_ARGS__); \
  }

#define VICOMH__(arg) vvi::detail::compare(lhs.arg, rhs.arg)

#define VVI_ICOM_EXPAND__(...) MKFNS(VVI_ICOM_EXPAND__,##__VA_ARGS__)
#define VVI_ICOM_EXPAND__0(void) true
#define VVI_ICOM_EXPAND__1(arg1, void) \
  VICOMH__(arg1)
#define VVI_ICOM_EXPAND__2(arg1, arg2, void) \
  VICOMH__(arg1) && VICOMH__(arg2)
#define VVI_ICOM_EXPAND__3(arg1, arg2, arg3, void) \
  VICOMH__(arg1) && VICOMH__(arg2) && VICOMH__(arg3)
#define VVI_ICOM_EXPAND__4(arg1, arg2, arg3, arg4, void) \
  VICOMH__(arg1) && VICOMH__(arg2) && VICOMH__(arg3) && \
  VICOMH__(arg4) 
#define VVI_ICOM_EXPAND__5(arg1, arg2, arg3, arg4, arg5, void) \
  VICOMH__(arg1) && VICOMH__(arg2) && VICOMH__(arg3) && \
  VICOMH__(arg4) && VICOMH__(arg5)
#define VVI_ICOM_EXPAND__6(arg1, arg2, arg3, arg4, arg5, arg6, void) \
  VICOMH__(arg1) && VICOMH__(arg2) && VICOMH__(arg3) &&\
  VICOMH__(arg4) && VICOMH__(arg5) && VICOMH__(arg6)
#define VVI_ICOM_EXPAND__7(arg1, arg2, arg3, arg4, arg5, arg6, arg7, void) \
  VICOMH__(arg1) && VICOMH__(arg2) && VICOMH__(arg3) && \
  VICOMH__(arg4) && VICOMH__(arg5) && VICOMH__(arg6) && \
  VICOMH__(arg7) 
#define VVI_ICOM_EXPAND__8(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, void) \
  VICOMH__(arg1) && VICOMH__(arg2) && VICOMH__(arg3) && \
  VICOMH__(arg4) && VICOMH__(arg5) && VICOMH__(arg6) && \
  VICOMH__(arg7) && VICOMH__(arg8)
#define VVI_ICOM_EXPAND__9(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, void) \
  VICOMH__(arg1) && VICOMH__(arg2) && VICOMH__(arg3) && \
  VICOMH__(arg4) && VICOMH__(arg5) && VICOMH__(arg6) && \
  VICOMH__(arg7) && VICOMH__(arg8) && VICOMH__(arg9)


// ----------------------------------------------------------------------------
// VALUE class declaration macros
// ----------------------------------------------------------------------------

#define VVI_VALUE_TOVIEW_OP__() \
  ViewType toView() { return ViewType(*this); }


#define VALUE_INTERFACE_METACLASS(view_t) \
public: \
  using ViewType = UNPACK view_t; \
  using ValueType = typename ViewType::ValueType; \
  using ImplType = typename ViewType::ImplType; \
  using ViewInterfaceType = typename ViewType::ViewInterfaceType; \
protected: \
  using Base = ::vvi::detail::ValueInterface<ViewType>; \
  using SmartPtrType = typename ViewType::SmartPtrType; \
public: \
  VVI_VALUE_TOVIEW_OP__() \
private:


#define VALUE_METACLASS_DECL_BEGIN(value_t, view_t) \
UNPACK value_t : public ::vvi::detail::ValueInterface<UNPACK view_t> { \
  VALUE_INTERFACE_METACLASS((UNPACK view_t)) \


#define VALUE_METACLASS_DECL_END() \
}



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

#define VVI_VIEW_DEF_CTOR(view_t) \
  view_t() = default;


#define VIEW_INTERFACE_METACLASS(value_t, view_t, impl_t) \
public: \
  using ImplType = UNPACK impl_t; \
  using ViewType = UNPACK view_t; \
  using ValueType = UNPACK value_t; \
  VVI_VIEW_DEF_CTOR(view_t) \
protected: \
  using Base = ::vvi::detail::ViewInterface<UNPACK impl_t>; \
  using ViewInterfaceType = Base; \
  using SmartPtrType = typename Base::SmartPtrType; \
  friend class UNPACK value_t; \
  VVI_VIEW_DEFINE_ALLOC_CTOR__(view_t) \
private:

#define VIEW_METACLASS_DECL_BEGIN(value_t, view_t, impl_t) \
UNPACK view_t : public ::vvi::detail::ViewInterface< UNPACK impl_t> { \
  VIEW_INTERFACE_METACLASS((UNPACK value_t), (UNPACK view_t), (UNPACK impl_t))

#define VIEW_METACLASS_DECL_END() \
}


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

