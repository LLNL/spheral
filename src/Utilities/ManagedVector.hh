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

template<typename U>
bool compare(ManagedVector<U> const&, ManagedVector<U> const&);

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
  {printf("MV Copy w/ Value Semantics.\n");
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

  template<typename U>
  friend bool compare(ManagedVector const& lhs, ManagedVector const& rhs);

};

template<typename U>
inline
bool compare(ManagedVector<U> const& lhs, ManagedVector<U> const& rhs)
{
  if (lhs.size() != rhs.size()) return false;
  for (size_t i = 0; i < lhs.size(); i++) {
    if (lhs[i] != rhs[i]) { 
      return false;
    }
  }
  return true;
}


template<typename U>
inline 
ManagedVector<U> deepCopy(ManagedVector<U> const& array)
{
  ManagedVector<U> copy(array.size());
  for (size_t i = 0; i < array.size(); i++) new (&copy[i]) U(array[i]);
  return copy;
}

} //  namespace Spheral

//#include "SphArrayInline.hh"

#else

//// Forward declare the SphArrayIterator class.
//namespace Spheral {
//  template<typename sph_array_t> class SphArrayIterator ;
//} //  namespace Spheral

#endif //  __Spheral_lvarray_hh__


