#ifndef __Spheral_lvarray_hh__
#define __Spheral_lvarray_hh__

#include "config.hh"
//#include "LvArray/Array.hpp"
//#include "LvArray/ChaiBuffer.hpp"
#include "chai/ManagedArray.hpp"

#include <cstdint>

constexpr uint32_t pow2_ceil(uint32_t v) {
    v--;
    v |= v >> 1;
    v |= v >> 2;
    v |= v >> 4;
    v |= v >> 8;
    v |= v >> 16;
    v++;

    return v;
}

namespace Spheral {

//#define SPHERAL_CALLBACK_ENABLED
#define STD_OUT_LOG

#ifdef STD_OUT_LOG
#define SPHERAL_LOG(A, B) std::cout << #A << " : " << B << std::endl;
#else
#define SPHERAL_LOG(A, B) UMPIRE_LOG(A, B)
#endif



//#define MV_VALUE_SEMANTICS

template<typename DataType>
class ManagedVector;

template<typename U>
ManagedVector<U> deepCopy(ManagedVector<U> const& array);

template<typename DataType>
class ManagedVector:
  private chai::ManagedArray<DataType>{
  using MA = chai::ManagedArray<DataType>;

public:

  inline static constexpr size_t initial_capacity = 8u;

  using MA::setUserCallback;

  using iterator = DataType*;
  using const_iterator = const DataType*;

  iterator begin() { return MA::begin(); }
  const_iterator begin() const { return MA::begin(); }

  iterator end() { return begin() + m_size; }
  const_iterator end() const { return begin() + m_size; }

  // ---------------------
  // Constructors
  // ---------------------
  SPHERAL_HOST_DEVICE ManagedVector() :
    MA()
  {
#if !defined(SPHERAL_GPU_ACTIVE) 
    setCallback();

#endif // SPHERAL_GPU_ACTIVE
  }

  SPHERAL_HOST_DEVICE ManagedVector(size_t elems) : 
    MA(),
    m_size(elems) 
  {
#if !defined(SPHERAL_GPU_ACTIVE) 
    MA::allocate(m_size < initial_capacity ? initial_capacity: pow2_ceil(m_size), chai::CPU, getCallback());
    for (size_t i = 0; i < m_size; i++) new (&MA::operator[](i)) DataType(); 
    MA::registerTouch(chai::CPU);
#endif // SPHERAL_GPU_ACTIVE
  }

  SPHERAL_HOST ManagedVector(size_t elems, DataType identity) :
    MA(),
    m_size(elems) 
  {
#if !defined(SPHERAL_GPU_ACTIVE) 
    MA::allocate(m_size < initial_capacity ? initial_capacity: pow2_ceil(m_size), chai::CPU, getCallback());
    for (size_t i = 0; i < m_size; i++) new (&MA::operator[](i)) DataType(identity);
    MA::registerTouch(chai::CPU);
#endif // SPHERAL_GPU_ACTIVE
  }

#ifdef MV_VALUE_SEMANTICS
  // ---------------------
  // Destructor
  // ---------------------
  SPHERAL_HOST ~ManagedVector() 
  {
    MA::free();
  }
#endif
  using MA::move;
  using MA::free;
  
  // ---------------------
  // Copy Constructor
  // ---------------------
#ifdef MV_VALUE_SEMANTICS
  SPHERAL_HOST_DEVICE constexpr inline ManagedVector(ManagedVector const& rhs) noexcept : 
    ManagedVector(rhs.m_size)
  {
    for (size_t i = 0; i < m_size; i++) new (&MA::operator[](i)) DataType(rhs[i]);
  }
#else
  SPHERAL_HOST_DEVICE constexpr inline ManagedVector(ManagedVector const& rhs) noexcept : MA(rhs), m_size(rhs.m_size) {}
#endif

  // ---------------------
  // Assignment
  // ---------------------
#ifdef MV_VALUE_SEMANTICS
  SPHERAL_HOST_DEVICE ManagedVector<DataType>& operator=(ManagedVector const& rhs) { 
    if (capacity() != rhs.capacity()) MA::reallocate(rhs.capacity());
    m_size = rhs.m_size;
    for (size_t i = 0; i < m_size; i++) new (&MA::operator[](i)) DataType(rhs[i]);
    return *this; 
  }
#else
  SPHERAL_HOST_DEVICE ManagedVector<DataType>& operator=(ManagedVector const& rhs) {
    MA::operator=(rhs);
    m_size = rhs.m_size;
    return *this; 
  }
#endif

  // ---------------------
  // Equivalence
  // ---------------------
#ifdef MV_VALUE_SEMANTICS
  SPHERAL_HOST bool operator==(ManagedVector const& rhs) const {
    if (m_size != rhs.m_size) return false;
    for (size_t i = 0; i < m_size; i++) {
      if (MA::operator[](i) != rhs[i]) { 
        return false;
      }
    }
    return true;
  }
#else
  SPHERAL_HOST_DEVICE bool operator==(ManagedVector const& rhs) const {
    if (m_size != rhs.m_size) return false;
    return MA::operator==(rhs);
  }
#endif
  SPHERAL_HOST_DEVICE bool operator!=(ManagedVector const& rhs) const {
    return !(*this == rhs);
  }

  SPHERAL_HOST void push_back(const DataType& value) {
    if (capacity() == 0) MA::allocate(initial_capacity, chai::CPU, getCallback());
    if (m_size >= capacity()) MA::reallocate(pow2_ceil(m_size + 1));
    new(&MA::operator[](m_size)) DataType(value);
    m_size++;
  }

  SPHERAL_HOST void push_back(DataType&& value) {
    if (capacity() == 0) MA::allocate(initial_capacity, chai::CPU, getCallback());
    if (m_size >= capacity()) MA::reallocate(pow2_ceil(m_size + 1));
    //MA::operator[](m_size) = std::move(value);
    new(&MA::operator[](m_size)) DataType(value);
    m_size++;
  }
  template<typename... Args>
  SPHERAL_HOST
  DataType& emplace_back(Args&&... args) {
    if (capacity() == 0) MA::allocate(initial_capacity, chai::CPU, getCallback());
    if (m_size >= capacity()) MA::reallocate(pow2_ceil(m_size + 1));

    new(&MA::data()[m_size]) DataType(std::forward<Args>(args)...);
    return MA::data()[m_size++];
  }

  SPHERAL_HOST
  void reserve(size_t c) {
    if (capacity() == 0) MA::allocate(c < initial_capacity ? initial_capacity: pow2_ceil(c), chai::CPU, getCallback());
    if (c >= capacity()) MA::reallocate(pow2_ceil(c));
  }

  SPHERAL_HOST
  void resize(size_t size) {
    const size_t old_size = m_size;

    if (old_size != size){
      if (old_size < size) {
        if (capacity() == 0) MA::allocate(size < initial_capacity ? initial_capacity: pow2_ceil(size), chai::CPU, getCallback());
        else if (capacity() < size) MA::reallocate(pow2_ceil(size));
        for (size_t i = old_size; i < size; i++) new(&MA::data(chai::CPU, false)[i]) DataType();
      }
      if (old_size > size) {
        destroy(begin() + old_size, begin() + size);
      }
      m_size = size;
      MA::registerTouch(chai::CPU);
    }
  }

  SPHERAL_HOST
  void insert(iterator pos, DataType const& value) {
    auto delta = std::distance(begin(), pos);
    if (m_size == 0) {
      push_back(value);
    }else {
      resize(m_size + 1);
      for (iterator it = end() - 1; it > begin() + delta; it--) {
        *it = std::move(*(it - 1));
        //*it = (*(it - 1));
      } 
      new(begin() + delta) DataType(value);
    }
  }

  SPHERAL_HOST
  void clear() {
    destroy(begin(), end());
    m_size = 0;
  }

  SPHERAL_HOST
  void erase(iterator pos) {
    for (iterator it = pos; it < end(); it++) {
      *it = std::move(*(it + 1));
    }
    m_size--;
  }

  SPHERAL_HOST_DEVICE size_t capacity() const {return MA::size();}
  SPHERAL_HOST_DEVICE size_t size() const {return m_size;}

  SPHERAL_HOST_DEVICE DataType& operator[](size_t idx) {return MA::data()[idx]; }
  SPHERAL_HOST_DEVICE DataType& operator[](size_t idx) const {return MA::data()[idx]; }


  // *******************************************************
  // Required to Allow ManagedVector to be properly CHAICopyable
  SPHERAL_HOST_DEVICE ManagedVector<DataType>& operator=(std::nullptr_t) { MA::operator=(nullptr); return *this; }
  SPHERAL_HOST_DEVICE void shallowCopy(const ManagedVector& other) {
    m_size=other.m_size;
    MA::shallowCopy(other);
  }
  // *******************************************************
  
  //SPHERAL_HOST operator std::vector<DataType>() const { return std::vector<DataType>(begin(), end()); }

  template< typename U=ManagedVector< DataType > >
  SPHERAL_HOST
  void setCallback() {
    MA::setUserCallback( getCallback() );
  }

  template< typename U=ManagedVector< DataType > >
  SPHERAL_HOST
  auto getCallback() {
//#ifdef SPHERAL_CALLBACK_ENABLED
//    std::string const typeString = LvArray::system::demangleType< U >();
//    return [typeString] (const chai::PointerRecord* record, chai::Action action, chai::ExecutionSpace exec) {
//        std::string const size = LvArray::system::calculateSize(record->m_size);
//        std::string const paddedSize = std::string( 9 - size.size(), ' ' ) + size;
//        char const * const spaceStr = ( exec == chai::CPU ) ? "HOST  " : "DEVICE";
//
//        if (action == chai::Action::ACTION_MOVE){
//          SPHERAL_LOG(Info, "Moved " << paddedSize << " to the " << spaceStr << ": " << typeString << " @ " <<  record->m_pointers[exec] )
//        }
//        if (action == chai::Action::ACTION_ALLOC){
//          SPHERAL_LOG(Info, "Allocated on " << spaceStr << " " << paddedSize << " : " << typeString << " @ " <<  record->m_pointers[exec] )
//        }
//        if (action == chai::Action::ACTION_FREE){
//          SPHERAL_LOG(Info, "Deallocated " << paddedSize << " : " << typeString << " @ " <<  record->m_pointers[exec] )
//        }
//      };
//#else
    return [](const chai::PointerRecord* , chai::Action , chai::ExecutionSpace ) {};
//#endif
  }

private:
  size_t m_size = 0;

  SPHERAL_HOST void destroy(iterator first, iterator last) {
    if ( !std::is_trivially_destructible< DataType >::value ) {
      for (iterator it = first; it < last; it++) {
        *it = DataType(); 
      }
    }
  }

  ManagedVector(MA const& managed_array) : MA(managed_array), m_size(managed_array.size()) {}

  template<typename U>
  friend ManagedVector deepCopy(ManagedVector const& array);

  SPHERAL_HOST_DEVICE
  friend bool compare(ManagedVector const& lhs, ManagedVector const& rhs)
  {
    if (lhs.m_size != rhs.m_size) return false;
    for (size_t i = 0; i < lhs.m_size; i++) {
      if (lhs[i] != rhs[i]) { 
        return false;
      }
    }
    return true;
  }

};


template<typename U>
inline 
ManagedVector<U> deepCopy(ManagedVector<U> const& array)
{
  ManagedVector<U> copy(array.size());
  for (size_t i = 0; i < array.size(); i++) new (&copy[i]) U(array[i]);
  return copy;
}

template<typename T>
class ManagedSmartPtr;

template<typename U>
ManagedSmartPtr<U> deepCopy(ManagedSmartPtr<U> const& rhs);

template<typename T>
ManagedSmartPtr<T> make_ManagedSmartPtr(T* host_ptr);

template<typename T, typename... Args>
ManagedSmartPtr<T> make_ManagedSmartPtr(Args... args);

template<typename T>
class ManagedSmartPtr : public chai::CHAICopyable
{
public:
  using element_type = T;
struct PrivateConstruct {};

  SPHERAL_HOST_DEVICE ManagedSmartPtr() {
#if !defined(SPHERAL_GPU_ACTIVE) 
    //m_ref_count = new size_t(0);
#endif // SPHERAL_GPU_ACTIVE
  }

  template<typename... Args>
  SPHERAL_HOST ManagedSmartPtr(PrivateConstruct, Args... args) {
#if !defined(SPHERAL_GPU_ACTIVE) 
    m_ptr = chai::ManagedArray<T>();
    m_ptr.allocate(1, chai::CPU, getCallback());
    m_ptr[0] = T(args...);
    m_ptr.registerTouch(chai::CPU);
    m_ref_count = new size_t(1);

#endif // SPHERAL_GPU_ACTIVE
  }

  SPHERAL_HOST ManagedSmartPtr(PrivateConstruct, T* host_ptr) {
#if !defined(SPHERAL_GPU_ACTIVE) 
    m_ptr = chai::makeManagedArray(host_ptr, 1, chai::CPU, true);
    m_ptr.setUserCallback(getCallback());
    m_ptr.registerTouch(chai::CPU);
    m_ref_count = new size_t(1);
#endif // SPHERAL_GPU_ACTIVE
  }


public:
  SPHERAL_HOST void registerTouch(chai::ExecutionSpace space) { m_ptr.registerTouch(space); }

  SPHERAL_HOST_DEVICE ManagedSmartPtr& operator=(ManagedSmartPtr const& rhs) {
    if (this != &rhs) {
      if (m_ptr != rhs.m_ptr) discontinue_ownership();
      m_ptr = rhs.m_ptr;
      m_ref_count = rhs.m_ref_count;
      increment_ref_count();
    }
    return *this;
  }

  SPHERAL_HOST_DEVICE ManagedSmartPtr(ManagedSmartPtr const& rhs) : m_ptr(rhs.m_ptr), m_ref_count(rhs.m_ref_count) {
    increment_ref_count();
  }

  SPHERAL_HOST void move(chai::ExecutionSpace space, bool touch = true) const { 
    m_ptr[0].move(space, touch);
  }

  SPHERAL_HOST_DEVICE T* get() const { return (m_ref_count) ? (m_ptr.data()) : nullptr; }
  SPHERAL_HOST_DEVICE T* operator->() { return get(); }
  SPHERAL_HOST_DEVICE T* operator->() const { return get(); }
  SPHERAL_HOST_DEVICE T& operator*() { return *get(); }
  SPHERAL_HOST_DEVICE T& operator*() const { return *get(); }

  SPHERAL_HOST_DEVICE ~ManagedSmartPtr() {
    discontinue_ownership();
  }

  SPHERAL_HOST_DEVICE ManagedSmartPtr& operator=(std::nullptr_t) { m_ref_count=nullptr; m_ptr=nullptr; return *this; }
  SPHERAL_HOST_DEVICE void shallowCopy(ManagedSmartPtr const& rhs) {
    *this = rhs;
  }

  template< typename U=ManagedSmartPtr< T > >
  SPHERAL_HOST
  auto getCallback() {
//#ifdef SPHERAL_CALLBACK_ENABLED
//    std::string const typeString = LvArray::system::demangleType< U >();
//    return [typeString] (const chai::PointerRecord* record, chai::Action action, chai::ExecutionSpace exec) {
//        std::string const size = LvArray::system::calculateSize(record->m_size);
//        std::string const paddedSize = std::string( 9 - size.size(), ' ' ) + size;
//        char const * const spaceStr = ( exec == chai::CPU ) ? "HOST  " : "DEVICE";
//
//        if (action == chai::Action::ACTION_MOVE){
//          SPHERAL_LOG(Info, "Moved " << paddedSize << " to the " << spaceStr << ": " << typeString << " @ " <<  record->m_pointers[exec] )
//        }
//        if (action == chai::Action::ACTION_ALLOC){
//          SPHERAL_LOG(Info, "Allocated on " << spaceStr << " " << paddedSize << " : " << typeString << " @ " <<  record->m_pointers[exec] )
//        }
//        if (action == chai::Action::ACTION_FREE){
//          SPHERAL_LOG(Info, "Deallocated " << paddedSize << " : " << typeString << " @ " <<  record->m_pointers[exec] )
//        }
//      };
//#else
    return [](const chai::PointerRecord* , chai::Action , chai::ExecutionSpace ) {};
//#endif
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
        m_ptr[0].free();
        m_ptr.free();
        delete m_ref_count;
        m_ref_count = nullptr;
      }
    }
#endif // SPHERAL_GPU_ACTIVE
  }


  SPHERAL_HOST_DEVICE
  friend bool compare(ManagedSmartPtr const& lhs, ManagedSmartPtr const& rhs)
  {
    // TODO : not safe
    return compare(lhs.m_ptr[0], rhs.m_ptr[0]);
  }

  chai::ManagedArray<T> m_ptr;
  size_t* m_ref_count = nullptr;

  template<typename U>
  friend ManagedSmartPtr deepCopy(ManagedSmartPtr const& rhs);

  template<typename U>
  friend ManagedSmartPtr make_ManagedSmartPtr(U* host_ptr);

  template<typename U, typename... Args>
  friend ManagedSmartPtr make_ManagedSmartPtr(Args... args);

};


template<typename U>
ManagedSmartPtr<U> make_ManagedSmartPtr(U* host_ptr)
  {
    ManagedSmartPtr<U> ptr = ManagedSmartPtr<U>(typename ManagedSmartPtr<U>::PrivateConstruct(), host_ptr);
    return ptr;
  }

template<typename U, typename... Args>
ManagedSmartPtr<U> make_ManagedSmartPtr(Args... args)
  {
    ManagedSmartPtr<U> ptr = ManagedSmartPtr<U>(typename ManagedSmartPtr<U>::PrivateConstruct(), args...);
    return ptr;
  }

template<typename U>
ManagedSmartPtr<U> deepCopy(ManagedSmartPtr<U> const& rhs)
{
  // TODO : not safe
  ManagedSmartPtr<U> ptr = make_ManagedSmartPtr<U>(deepCopy(*rhs));
  return ptr;
  //return ManagedSmartPtr<U>(typename ManagedSmartPtr<U>::PrivateConstruct(), deepCopy(*rhs));
}


template<typename T>
class MVSmartRef : public ManagedSmartPtr<ManagedVector<T>>
{
  using Base = ManagedSmartPtr<ManagedVector<T>>;
  SPHERAL_HOST_DEVICE ManagedVector<T> & mv() { return Base::m_ptr[0]; }
  SPHERAL_HOST_DEVICE ManagedVector<T> const& mv() const { return Base::m_ptr[0]; }

public:

  using Base::operator->;
  using Base::operator*;
  using Base::get;
  using Base::move;

  using MV = Spheral::ManagedVector<T>;

  SPHERAL_HOST_DEVICE MVSmartRef() = default;

  template<typename... Args>
  MVSmartRef(Args... args) : Base(make_ManagedSmartPtr<MV>(args...)) {mv().setCallback();}
  
  using iterator = typename MV::iterator;
  using const_iterator = typename MV::const_iterator;

  iterator begin() { return mv().begin(); }
  const_iterator begin() const { return mv().begin(); }

  iterator end() { return begin() + size(); }
  const_iterator end() const { return begin() + size(); }

  SPHERAL_HOST_DEVICE T& operator[](size_t idx) {return mv()[idx]; }
  SPHERAL_HOST_DEVICE T& operator[](size_t idx) const {return mv()[idx]; }

  SPHERAL_HOST_DEVICE size_t size() const { return mv().size(); }

  SPHERAL_HOST void resize(size_t sz) {
    move(chai::CPU);
    mv().resize(sz);
    Base::m_ptr.registerTouch(chai::CPU);
  }

  SPHERAL_HOST
  void insert(iterator pos, T const& value) {
    move(chai::CPU);
    mv().insert(pos, value);
    Base::m_ptr.registerTouch(chai::CPU);
  }

  SPHERAL_HOST void push_back(T const& value) {
    move(chai::CPU);
    mv().push_back(value);
    Base::m_ptr.registerTouch(chai::CPU);
  }
  
  SPHERAL_HOST void push_back(T&& value) {
    move(chai::CPU);
    mv().push_back(value);
    Base::m_ptr.registerTouch(chai::CPU);
  }
  
  SPHERAL_HOST
  void reserve(size_t c) {
    move(chai::CPU);
    mv().reserve(c);
    Base::m_ptr.registerTouch(chai::CPU);
  }
  
  SPHERAL_HOST
  void clear() {
    move(chai::CPU);
    mv().clear();
    Base::m_ptr.registerTouch(chai::CPU);
  }

  SPHERAL_HOST
  void erase(iterator pos) {
    move(chai::CPU);
    mv().erase(pos);
    Base::m_ptr.registerTouch(chai::CPU);
  }

  SPHERAL_HOST_DEVICE MVSmartRef& operator=(std::nullptr_t) { Base::operator=(nullptr); return *this; }
  SPHERAL_HOST_DEVICE void shallowCopy(MVSmartRef const& rhs) {
    Base::shallowCopy(rhs);
  }

  SPHERAL_HOST_DEVICE bool operator==(MVSmartRef const& rhs) const { return (mv() == rhs.mv()); }
  SPHERAL_HOST_DEVICE bool operator!=(MVSmartRef const& rhs) const { return (mv() != rhs.mv()); }

private:

  friend MVSmartRef deepCopy(MVSmartRef const& rhs)
  {
    // TODO : not safe
    return MVSmartRef(deepCopy(rhs.m_ptr[0]));
  }

  SPHERAL_HOST_DEVICE
  friend bool compare(MVSmartRef const& lhs, MVSmartRef const& rhs)
  {
    // TODO : not safe
    return compare(lhs.mv(), rhs.mv());
  }
};




//template<typename DataType>
//using SphArray = LvArray::Array<DataType, 1, camp::idx_seq<0>, std::ptrdiff_t, LvArray::ChaiBuffer>;
//
//template<typename DataType>
//using SphArrayView = LvArray::ArrayView<DataType, 1, 0, std::ptrdiff_t, LvArray::ChaiBuffer>;

//template<typename sph_array_t>
//class SphArrayIterator {
//public:
//  using iterator_category = std::random_access_iterator_tag;
//  using value_type = typename sph_array_t::ValueType;
//  using difference_type = std::ptrdiff_t;
//  using pointer = value_type*;
//  using reference = value_type&;
//  
//  SphArrayIterator(pointer ptr);
//
//  SphArrayIterator& operator++();
//  SphArrayIterator operator++(int);
//
//  SphArrayIterator& operator--();
//  SphArrayIterator operator--(int);
//
//  SphArrayIterator operator+(const difference_type& index) const;
//  SphArrayIterator operator-(const difference_type& index) const;
//
//  difference_type operator-(const SphArrayIterator& it);
//
//  reference operator[](int index);
//  pointer operator->();
//
//  reference operator*() const;
//
//  bool operator==(const SphArrayIterator& rhs) const;
//  bool operator!=(const SphArrayIterator& rhs) const;
//  bool operator<(const SphArrayIterator& rhs) const;
//  bool operator<=(const SphArrayIterator& rhs) const;
//  bool operator>(const SphArrayIterator& rhs) const;
//  bool operator>=(const SphArrayIterator& rhs) const;
//
//  operator SphArrayIterator<typename sph_array_t::ViewTypeConst>() const;
//
//private:
//  pointer mPtr;
//};
//
//template<typename sph_array_t>
//class SphArrayFieldIterator {
//public:
//  using iterator_category = std::input_iterator_tag;
//  using value_type = typename sph_array_t::ValueType;
//  using difference_type = std::ptrdiff_t;
//  using pointer = value_type*;
//  using reference = value_type&;
//
//  using field_type = typename value_type::FieldType;
//  using field_pointer = field_type*;
//  using field_reference = field_type&;
//
//  SphArrayFieldIterator(pointer ptr);
//
//  field_reference operator*() const;
//  field_pointer operator->();
//
//  SphArrayFieldIterator& operator++();
//  SphArrayFieldIterator operator++(int);
//
//  bool operator==(const SphArrayFieldIterator& rhs) const;
//
//private:
//  pointer mPtr;
//};


} //  namespace Spheral

//#include "SphArrayInline.hh"

#else

//// Forward declare the SphArrayIterator class.
//namespace Spheral {
//  template<typename sph_array_t> class SphArrayIterator ;
//} //  namespace Spheral

#endif //  __Spheral_lvarray_hh__


