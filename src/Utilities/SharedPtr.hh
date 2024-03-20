#ifndef __Spheral_SharedPtr_hh__
#define __Spheral_SharedPtr_hh__

#include "config.hh"
#include "chai/ManagedArray.hpp"
#include <memory>
#include <type_traits>

namespace Spheral {


// std::shared_ptr supports a variety of features that I have chosen not 
// to support in this implementation:
//     1. Array capture: std::shared_ptr since c++17 can manage a bounded or unbounded array. this
//        implementation only manages single objects.
//     2. Thread safety: std::shared_ptr maintains thread safety through the use of mutex symantics 
//        and atomic types.
//     3. Support for enable_shared_from_this. Which is a way to have a facotry class only create
//        objects of itself as a shared_ptr.
//     4. Conversion form std::weak_ptr and std::unique_ptr.
//     5. EBO (Empty Base Optimization) / EBCO (Class): std::shared_ptr takes advatage of EBO in 
//        it's implementation of the deleter counter. We are doing a much simpler implementation
//        here.
//     6. Optional Allocator: since we are not supporting arrays this is mostly useless to us. It 
//        also means we will not support allocate_shared().

// _Sp_counted_base is the base type for shared_ptr's reference counting funcitonality.
// Since shared_ptr is responsible for owning and deleteing the underlying memory of an 
// object the reference count object will actually handle this, shared_ptr does not deallocate 
// the object itself.

// std::shared_ptr holds a pointer to a base reference counter so that it can use type-erasure
// on the specific counter type a user defines. The counter type is specified by the Ctor arguments
// a user provides. The counter type also determines the behavior of the dispose and destroy 
// functionality of the owned pointer. For this it is importnt to us to be able to pass a custom
// deleter to our Spheral::shared_ptr.


// ****************************************************************************
// _Sp_counted_base
// 
// The std derives this from mutex like classes that are template selectable on
// the locking behavior/ policy. We are not implementing that here...
class _Sp_counted_base {
public:
  _Sp_counted_base() noexcept : _M_use_count(1) {}

  virtual ~_Sp_counted_base() noexcept {}
  
  virtual void _M_dispose() noexcept = 0;

  virtual void _M_destroy() noexcept { delete this; }

  void _M_add_ref_lock() {
    assert(_M_add_ref_lock_nothrow()); // would usually throw std::bad_weak_ptr
  }
  
  bool _M_add_ref_lock_nothrow() noexcept {
    if (_M_use_count == 0) return false;
    ++_M_use_count;
    return true;
  }

  void _M_add_ref_copy() noexcept { ++_M_use_count; }

  void _M_release() noexcept {
    if (--_M_use_count == 0) {
      _M_dispose();
      _M_destroy(); // Would check weak_count before doing this...
    }
  }

  long _M_get_use_count() const noexcept { return _M_use_count; } 

private:
  _Sp_counted_base(_Sp_counted_base const&) = delete;
  _Sp_counted_base& operator=(_Sp_counted_base const&) = delete;

  long _M_use_count = 0; // Would be some atomic integer type.
  // Not counting weak_count, as we aren't implementing optimizations for make_shared yet.
  // https://stackoverflow.com/questions/49585818/why-does-shared-ptr-needs-to-hold-reference-counting-for-weak-ptr
};


// ****************************************************************************
// _Sp_counted_ptr
//
// This is the "default" counter object. It holds the owend pointer and will call a simple delete 
// when the use_count == 0, not much to explain here.
template<typename _Ptr>
class _Sp_counted_ptr final : public _Sp_counted_base {
public:

  _Sp_counted_ptr(_Ptr p) noexcept : _M_ptr(p) {}

  virtual void _M_dispose() noexcept { delete _M_ptr; }

  virtual void _M_destroy() noexcept { delete this; }

  _Sp_counted_ptr(_Sp_counted_ptr const&) = delete;
  _Sp_counted_ptr& operator=(_Sp_counted_ptr const&) = delete;

private:

  _Ptr _M_ptr;
};



// ****************************************************************************
// _Sp_counted_deleter
//
template<typename _Ptr, typename _Deleter>
class _Sp_counted_deleter final : public _Sp_counted_base {

  // counted_deleter typically uses EBO helpers here to disgues the extra pointer in the object
  // size, it also accepts allocators which we aren't doing in this implementation.
  class _Impl {
    public:
      _Impl(_Ptr p, _Deleter d) : _M_ptr(p), _M_deleter(std::move(d)) {}

      _Deleter& _M_del() noexcept { return _M_deleter; }

      _Ptr _M_ptr;
      _Deleter _M_deleter;
  };

public:
  _Sp_counted_deleter(_Ptr p, _Deleter d) : _M_impl(p, std::move(d)) {}
  ~_Sp_counted_deleter() noexcept {}

  virtual void _M_dispose() noexcept { _M_impl._M_del()(_M_impl._M_ptr); }
  virtual void _M_destroy() noexcept { this->~_Sp_counted_deleter(); }

private:
  _Impl _M_impl;
};


// ****************************************************************************
// __shared_count
//
// __shared_count is an interface class used to define the which counter type the shared_ptr
// needs to construct. This is done through Ctor argument overloads. In this case when 
// a deleter argument is passed to the __shared_count Ctor it will create a _Sp_counted_deleter
//
// __shared_count is also responsible for calling _M_add_ref_copy() and _M_release() during 
// copy, assignemt and destruction of a __shared_count object. Removing that responsibility 
// from the higher level interfaces of shared_ptr.
class __shared_count {

public:
  constexpr __shared_count() noexcept : _M_pi(0) {}

  // TODO: Change this to use the new CTor explicit call not this heap allocation of new we have here.
  template<typename _Ptr>
  explicit __shared_count(_Ptr p) : _M_pi(new _Sp_counted_ptr<_Ptr>(p)) {}

  template<typename _Ptr, typename _Deleter>
  explicit __shared_count(_Ptr p, _Deleter d) : _M_pi(new _Sp_counted_deleter<_Ptr, _Deleter>(p, d)) {}

  ~__shared_count() noexcept { if(_M_pi) _M_pi->_M_release(); }

  __shared_count(__shared_count const& r) noexcept : _M_pi(r._M_pi) { if (_M_pi) _M_pi->_M_add_ref_copy(); }

  __shared_count& operator=(__shared_count const& r) noexcept {
    _Sp_counted_base* temp = r._M_pi;
    if (temp != _M_pi)
    {
      if (temp) temp->_M_add_ref_copy();
      if (_M_pi) _M_pi->_M_release();
      _M_pi = temp;
    }
    return *this;
  }

  void _M_swap(__shared_count& r) noexcept {
    _Sp_counted_base* temp = r._M_pi;
    r._M_pi = _M_pi;
    _M_pi = temp;
  }

  long _M_get_use_count() const noexcept { _M_pi ? _M_pi->_M_get_use_count() : 0; }

  friend inline bool
  operator==(__shared_count const& a, __shared_count const& b) noexcept { return a._M_pi == b._M_pi; }

  _Sp_counted_base* _M_pi;
};



// ****************************************************************************
// __shared_ptr
//
// The std separates shared_ptr from some of it's implementation in __shared_ptr.
// In this implementation I'm not certain how useful it is as std::__shared_ptr 
// does a lot more behind the scenes work for functionality Spheral::shared_ptr
// is not supporting. However it still helps to keep the interface class of 
// Spheral::shared_ptr very clean and readable.


// Type traits for SFINAE in __shared_ptr

// Trait to check if shared_ptr<T> can be constructed from Y*.
// Y* shall be convertible to T*.
template<typename _Tp, typename _Yp>
struct __sp_is_constructible : std::is_convertible<_Yp*, _Tp*>::type {};

// A pointer type Y* is said to be compatible with a pointer type T* when
// either Y* is convertible to T* or Y is U[N] and T is U cv [].
template<typename _Yp_ptr, typename _Tp_ptr>
struct __sp_compatible_with : std::false_type {};

template<typename _Yp, typename _Tp>
struct __sp_compatible_with<_Yp*, _Tp*> : std::is_convertible<_Yp*, _Tp*>::type {};



template<typename _Tp>
class __shared_ptr { 

public:
  using element_type = typename std::remove_extent<_Tp>::type;
  
private:
  // Constraint for taking ownership of a pointer of type _Yp*:
  template<typename _Yp>
  using _SafeConv = typename std::enable_if<__sp_is_constructible<_Tp, _Yp>::value>::type;

  // Constraint for construction from shared_ptr and weak_ptr:
  template<typename _Yp, typename _Res = void>
  using _Compatible = typename std::enable_if<__sp_compatible_with<_Yp*, _Tp*>::value, _Res>::type;

  // Constraint for assignment from shared_ptr and weak_ptr:
  template<typename _Yp>
  using _Assignable = _Compatible<_Yp, __shared_ptr&>;

public:

  constexpr __shared_ptr() noexcept : _M_ptr(0), _M_refcount() {}

  template<typename _Yp, typename = _SafeConv<_Yp>>
  explicit __shared_ptr(_Yp* p) : _M_ptr(p), _M_refcount(p) {}

  template<typename _Yp, typename _Deleter, typename = _SafeConv<_Yp>>
  __shared_ptr(_Yp* p, _Deleter d) : _M_ptr(p), _M_refcount(p, std::move(d)) {
    static_assert(std::is_invocable<_Deleter&, _Yp*&>::value, "deleter expression d(p) is well-formed");
  }

  __shared_ptr(const __shared_ptr&) noexcept = default;
  __shared_ptr& operator=(const __shared_ptr&) noexcept = default;
  ~__shared_ptr() = default;

  template<typename _Yp, typename = _Compatible<_Yp>>
	__shared_ptr(const __shared_ptr<_Yp>& __r) noexcept : _M_ptr(__r._M_ptr), _M_refcount(__r._M_refcount) {} 

  __shared_ptr(__shared_ptr&& __r) noexcept : _M_ptr(__r._M_ptr), _M_refcount() {
    _M_refcount._M_swap(__r._M_refcount);
    __r._M_ptr = nullptr;
  }

  template<typename _Yp, typename = _Compatible<_Yp>>
	__shared_ptr(__shared_ptr<_Yp>&& __r) noexcept : _M_ptr(__r._M_ptr), _M_refcount() {
	  _M_refcount._M_swap(__r._M_refcount);
	  __r._M_ptr = nullptr;
	}

  __shared_ptr& operator=(__shared_ptr&& __r) noexcept { __shared_ptr(std::move(__r)).swap(*this); return *this; }

  template<class _Yp>
	_Assignable<_Yp> operator=(__shared_ptr<_Yp>&& __r) noexcept {
	  __shared_ptr(std::move(__r)).swap(*this);
	  return *this;
	}

  void reset() noexcept { __shared_ptr().swap(*this); }

  template<typename _Yp>
	_SafeConv<_Yp> reset(_Yp* __p) // _Yp must be complete.
	{
	  // Catch self-reset errors.
	  assert(__p == nullptr || __p != _M_ptr);
	  __shared_ptr(__p).swap(*this);
	}

  template<typename _Yp, typename _Deleter>
	_SafeConv<_Yp> reset(_Yp* __p, _Deleter __d) { __shared_ptr(__p, std::move(__d)).swap(*this); }

  /// Return the stored pointer.
  element_type* get() const noexcept { return _M_ptr; }

  /// Return true if the stored pointer is not null.
  explicit operator bool() const noexcept { return _M_ptr != nullptr; }

  /// If *this owns a pointer, return the number of owners, otherwise zero.
  long use_count() const noexcept { return _M_refcount._M_get_use_count(); }

  /// Exchange both the owned pointer and the stored pointer.
  void swap(__shared_ptr<_Tp>& __other) noexcept
  {
    std::swap(_M_ptr, __other._M_ptr);
    _M_refcount._M_swap(__other._M_refcount);
  }

  element_type& operator*() const noexcept { assert(_M_get() != nullptr); return *_M_get(); }

  element_type* operator->() const noexcept { assert(_M_get() != nullptr); return _M_get(); }

private:
  element_type* _M_get() const noexcept { return static_cast<const __shared_ptr<_Tp>*>(this)->get(); }

  template<typename _Tp1>
  friend class __shared_ptr;

  element_type*	 _M_ptr;         // Contained pointer.
  __shared_count _M_refcount;    // Reference counter.
};



template<typename _Tp>
class shared_ptr : public __shared_ptr<_Tp> {

  template<typename... _Args>
	using _Constructible = typename std::enable_if< std::is_constructible<__shared_ptr<_Tp>, _Args...>::value >::type;

  template<typename _Arg>
	using _Assignable = typename std::enable_if< std::is_assignable<__shared_ptr<_Tp>&, _Arg>::value, shared_ptr& >::type;

public:

  // The type pointed to by the stored pointer, remove_extent_t<_Tp>
  using element_type = typename __shared_ptr<_Tp>::element_type;

  constexpr shared_ptr() noexcept : __shared_ptr<_Tp>() { }

  shared_ptr(const shared_ptr&) noexcept = default; ///< Copy constructor

  template<typename _Yp, typename = _Constructible<_Yp*>>
	explicit shared_ptr(_Yp* __p) : __shared_ptr<_Tp>(__p) { }

  template<typename _Yp, typename _Deleter, typename = _Constructible<_Yp*, _Deleter>>
	shared_ptr(_Yp* __p, _Deleter __d) : __shared_ptr<_Tp>(__p, std::move(__d)) { }

  template<typename _Yp, typename = _Constructible<const shared_ptr<_Yp>&>>
	shared_ptr(const shared_ptr<_Yp>& __r) noexcept : __shared_ptr<_Tp>(__r) { }

  shared_ptr(shared_ptr&& __r) noexcept : __shared_ptr<_Tp>(std::move(__r)) { }

  template<typename _Yp, typename = _Constructible<shared_ptr<_Yp>>>
	shared_ptr(shared_ptr<_Yp>&& __r) noexcept : __shared_ptr<_Tp>(std::move(__r)) { }


  shared_ptr& operator=(const shared_ptr&) noexcept = default;

  template<typename _Yp>
	_Assignable<const shared_ptr<_Yp>&> operator=(const shared_ptr<_Yp>& __r) noexcept
	{
	  this->__shared_ptr<_Tp>::operator=(__r);
	  return *this;
	}

  shared_ptr& operator=(shared_ptr&& __r) noexcept {
    this->__shared_ptr<_Tp>::operator=(std::move(__r));
    return *this;
  }

  template<class _Yp>
	_Assignable<shared_ptr<_Yp>> operator=(shared_ptr<_Yp>&& __r) noexcept {
	  this->__shared_ptr<_Tp>::operator=(std::move(__r));
	  return *this;
	}

}; // class shared_ptr


/// Equality operator for shared_ptr objects, compares the stored pointers
template<typename _Tp, typename _Up>
inline bool operator==(const shared_ptr<_Tp>& __a, const shared_ptr<_Up>& __b) noexcept { return __a.get() == __b.get(); }

template<typename _Tp, typename _Up>
inline bool operator!=(const shared_ptr<_Tp>& __a, const shared_ptr<_Up>& __b) noexcept { return __a.get() != __b.get(); }






template<typename T>
class managed_shared_ptr;

template<typename U>
managed_shared_ptr<U> deepCopy(managed_shared_ptr<U> const& rhs);

template<typename T, typename... Args>
managed_shared_ptr<T> make_managed_shared_ptr(Args... args);

template<typename T>
class managed_shared_ptr : public chai::CHAICopyable
{

  template<typename... Args>
    using Constructible = typename std::enable_if<
      std::is_constructible<managed_shared_ptr<T>, Args...>::value
    >::type;


public:
  using element_type = T;
struct PrivateConstruct {};

  SPHERAL_HOST_DEVICE managed_shared_ptr() {
#if !defined(SPHERAL_GPU_ACTIVE) 
    //m_ref_count = new size_t(0);
#endif // SPHERAL_GPU_ACTIVE
  }

  template<typename... Args>
  SPHERAL_HOST managed_shared_ptr(PrivateConstruct, Args... args) {
#if !defined(SPHERAL_GPU_ACTIVE) 
    m_ptr = chai::ManagedArray<T>();
    m_ptr.allocate(1, chai::CPU, getCallback());
    m_ptr[0] = T(args...);
    m_ptr.registerTouch(chai::CPU);
    m_ref_count = new size_t(1);
#endif // SPHERAL_GPU_ACTIVE
  }


  SPHERAL_HOST managed_shared_ptr(T* host_ptr) {
#if !defined(SPHERAL_GPU_ACTIVE) 
    m_ptr = chai::makeManagedArray(host_ptr, 1, chai::CPU, true);
    m_ptr.setUserCallback(getCallback());
    m_ptr.registerTouch(chai::CPU);
    m_ref_count = new size_t(1);
#endif // SPHERAL_GPU_ACTIVE
  }

  template<typename Y, typename = Constructible<Y*> >
  SPHERAL_HOST managed_shared_ptr(Y* host_ptr) : managed_shared_ptr(static_cast<T*>(host_ptr)) { }



public:
  SPHERAL_HOST void registerTouch(chai::ExecutionSpace space) { m_ptr.registerTouch(space); }

  SPHERAL_HOST_DEVICE managed_shared_ptr& operator=(managed_shared_ptr const& rhs) {
    if (this != &rhs) {
      if (m_ptr != rhs.m_ptr) discontinue_ownership();
      m_ptr = rhs.m_ptr;
      m_ref_count = rhs.m_ref_count;
      increment_ref_count();
    }
    return *this;
  }

  SPHERAL_HOST_DEVICE managed_shared_ptr(managed_shared_ptr const& rhs) : m_ptr(rhs.m_ptr), m_ref_count(rhs.m_ref_count) {
    increment_ref_count();
  }

  template<typename Y, typename = Constructible<Y*> >
	managed_shared_ptr(const managed_shared_ptr<Y>& rhs) noexcept : m_ptr(rhs.m_ptr), m_ref_count(rhs.m_ref_count) {
    increment_ref_count();
  }

  SPHERAL_HOST void move(chai::ExecutionSpace space, bool touch = true) const { 
    m_ptr[0].move(space, touch);
  }

  SPHERAL_HOST size_t use_count() const { return (m_ref_count) ? *m_ref_count : 0; }

  SPHERAL_HOST_DEVICE T* get() const { return (m_ref_count) ? (m_ptr.data()) : nullptr; }
  SPHERAL_HOST_DEVICE T* operator->() { return get(); }
  SPHERAL_HOST_DEVICE T* operator->() const { return get(); }
  SPHERAL_HOST_DEVICE T& operator*() { return *get(); }
  SPHERAL_HOST_DEVICE T& operator*() const { return *get(); }

  SPHERAL_HOST_DEVICE ~managed_shared_ptr() {
    discontinue_ownership();
  }

  SPHERAL_HOST_DEVICE managed_shared_ptr& operator=(std::nullptr_t) { m_ref_count=nullptr; m_ptr=nullptr; return *this; }
  SPHERAL_HOST_DEVICE void shallowCopy(managed_shared_ptr const& rhs) {
    *this = rhs;
  }

  template< typename U=managed_shared_ptr< T > >
  SPHERAL_HOST
  auto getCallback() {
#ifdef SPHERAL_CALLBACK_ENABLED
    std::string const typeString = LvArray::system::demangleType< U >();
    return [typeString] (const chai::PointerRecord* record, chai::Action action, chai::ExecutionSpace exec) {
        std::string const size = LvArray::system::calculateSize(record->m_size);
        std::string const paddedSize = std::string( 9 - size.size(), ' ' ) + size;
        char const * const spaceStr = ( exec == chai::CPU ) ? "HOST  " : "DEVICE";

        if (action == chai::Action::ACTION_MOVE){
          SPHERAL_LOG(Info, "Moved " << paddedSize << " to the " << spaceStr << ": " << typeString << " @ " <<  record->m_pointers[exec] )
        }
        if (action == chai::Action::ACTION_ALLOC){
          SPHERAL_LOG(Info, "Allocated on " << spaceStr << " " << paddedSize << " : " << typeString << " @ " <<  record->m_pointers[exec] )
        }
        if (action == chai::Action::ACTION_FREE){
          SPHERAL_LOG(Info, "Deallocated " << paddedSize << " : " << typeString << " @ " <<  record->m_pointers[exec] )
        }
      };
#else
    return [](const chai::PointerRecord* , chai::Action , chai::ExecutionSpace ) {};
#endif
  }
protected:

  SPHERAL_HOST_DEVICE void increment_ref_count() {
#if !defined(SPHERAL_GPU_ACTIVE) 
    if (m_ref_count != nullptr) (*m_ref_count)++;
#endif // SPHERAL_GPU_ACTIVE
  }

  SPHERAL_HOST_DEVICE void discontinue_ownership() {
#if !defined(SPHERAL_GPU_ACTIVE) 
    if (m_ref_count != nullptr){
      (*m_ref_count)--;
      if (*m_ref_count == 0)
      {
        //Should use a deleter here...
        //m_ptr[0].free();
        delete &m_ptr[0];
        //m_ptr.free();
        delete m_ref_count;
        m_ref_count = nullptr;
      }
    }
#endif // SPHERAL_GPU_ACTIVE
  }


  SPHERAL_HOST_DEVICE
  friend bool compare(managed_shared_ptr const& lhs, managed_shared_ptr const& rhs)
  {
    // TODO : not safe
    return compare(lhs.m_ptr[0], rhs.m_ptr[0]);
  }

public:
  chai::ManagedArray<T> m_ptr;
  size_t* m_ref_count = nullptr;

  template<typename U>
  friend managed_shared_ptr deepCopy(managed_shared_ptr const& rhs);

  template<typename U, typename... Args>
  friend managed_shared_ptr make_managed_shared_ptr(Args... args);
};

template<typename U, typename... Args>
managed_shared_ptr<U> make_managed_shared_ptr(Args... args)
  {
    managed_shared_ptr<U> ptr = managed_shared_ptr<U>(typename managed_shared_ptr<U>::PrivateConstruct(), args...);
    return ptr;
  }

template<typename U>
managed_shared_ptr<U> deepCopy(managed_shared_ptr<U> const& rhs)
{
  // TODO : not safe
  managed_shared_ptr<U> ptr = make_managed_shared_ptr<U>(deepCopy(*rhs));
  return ptr;
  //return managed_shared_ptr<U>(typename managed_shared_ptr<U>::PrivateConstruct(), deepCopy(*rhs));
}

} // namespace Spheral

#endif // __Spheral_SharedPtr_hh__
