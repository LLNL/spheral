#ifndef __Spheral_lvarray_hh__
#define __Spheral_lvarray_hh__

#include "LvArray/Array.hpp"
#include "LvArray/ChaiBuffer.hpp"
#include "chai/ManagedArray.hpp"

namespace Spheral {


template<typename DataType>
class ManagedVector:
  private chai::ManagedArray<DataType>{
  using MA = chai::ManagedArray<DataType>;

public:

  using iterator = DataType*;
  using const_iterator = const DataType*;

  iterator begin() { return MA::begin(); }
  const_iterator begin() const { return MA::begin(); }

  iterator end() { return begin() + m_size; }
  const_iterator end() const { return begin() + m_size; }

  // ---------------------
  // Constructors
  // ---------------------
  RAJA_HOST ManagedVector() : MA(6, chai::CPU) {}
  RAJA_HOST ManagedVector(size_t elems) : MA(elems, chai::CPU), m_size(elems) {
    for (size_t i = 0; i < m_size; i++) new (&MA::operator[](i)) DataType(); 
  }
  RAJA_HOST ManagedVector(size_t elems, DataType identity) : MA(elems, chai::CPU), m_size(elems) {
    for (size_t i = 0; i < m_size; i++) new (&MA::operator[](i)) DataType(identity);
  }
  
  // ---------------------
  // Copy Constructor
  // ---------------------
  RAJA_HOST_DEVICE constexpr inline ManagedVector(ManagedVector const& rhs) noexcept : MA(rhs.capacity()), m_size(rhs.m_size) {
    for (size_t i = 0; i < m_size; i++) new (&MA::operator[](i)) DataType(rhs[i]);
  }

  // ---------------------
  // Assignment
  // ---------------------
  RAJA_HOST_DEVICE ManagedVector<DataType>& operator=(ManagedVector const& rhs) { 
    if (capacity() != rhs.capacity()) MA::reallocate(rhs.capacity());
    m_size = rhs.m_size;
    for (size_t i = 0; i < m_size; i++) new (&MA::operator[](i)) DataType(rhs[i]);
    return *this; 
  }

  // ---------------------
  // Equivalence
  // ---------------------
  bool operator==(ManagedVector const& rhs) = delete;
  RAJA_HOST_DEVICE bool operator==(ManagedVector const& rhs) const {
    if (m_size != rhs.m_size) return false;
    for (size_t i = 0; i < m_size; i++) {
      if (MA::operator[](i) != rhs[i]) { 
        return false;
      }
    }
    return true;
  }

  RAJA_HOST void push_back(const DataType& value) {
    if (m_size >= capacity()) MA::reallocate(capacity() + (capacity() / 2));
    MA::data()[m_size] = value;
    m_size++;
  }
  RAJA_HOST void push_back(DataType&& value) {
    if (m_size >= capacity()) MA::reallocate(capacity() + (capacity() / 2));

    MA::data()[m_size] = std::move(value);
    m_size++;
  }
  template<typename... Args>
  RAJA_HOST
  DataType& emplace_back(Args&&... args) {
    if (m_size >= capacity()) MA::reallocate(capacity() + (capacity() / 2));

    new(&MA::data()[m_size]) DataType(std::forward<Args>(args)...);
    return MA::data()[m_size++];
  }

  RAJA_HOST
  void resize(size_t size) {
    const size_t old_size = m_size;

    if (old_size < size) {
      MA::reallocate(size);
      for (size_t i = old_size; i < size; i++) new(&MA::operator[](i)) DataType();
    }
    if (old_size > size) {
      destroy(begin() + old_size, begin() + size);
    }

    m_size = size;
  }

  RAJA_HOST
  void erase(iterator pos) {
    for (iterator it = pos; it < end(); it++) {
      *it = std::move(*(it + 1));
    }
    m_size--;
  }

  RAJA_HOST_DEVICE size_t capacity() const {return MA::m_elems;}
  RAJA_HOST_DEVICE size_t size() const {return m_size;}

  RAJA_HOST_DEVICE DataType& operator[](size_t idx) {return MA::data()[idx]; }
  RAJA_HOST_DEVICE DataType& operator[](size_t idx) const {return MA::data()[idx]; }


  // *******************************************************
  // Required to Allow ManagedVector to be properly CHAICopyable
  RAJA_HOST_DEVICE ManagedVector<DataType>& operator=(std::nullptr_t) { MA::operator=(nullptr); return *this; }
  RAJA_HOST_DEVICE void shallowCopy(const ManagedVector& other) {
    m_size=other.m_size;
    MA::shallowCopy(other);
  }
  // *******************************************************
  
  //RAJA_HOST operator std::vector<DataType>() const { return std::vector<DataType>(begin(), end()); }

private:
  size_t m_size = 0;

  RAJA_HOST void destroy(iterator first, iterator last) {
    if ( !std::is_trivially_destructible< DataType >::value ) {
      for (iterator it = first; it < last; it++) {
        *it = DataType(); 
      }
    }
  }

};

template<typename DataType>
using SphArray = LvArray::Array<DataType, 1, camp::idx_seq<0>, std::ptrdiff_t, LvArray::ChaiBuffer>;

template<typename DataType>
using SphArrayView = LvArray::ArrayView<DataType, 1, 0, std::ptrdiff_t, LvArray::ChaiBuffer>;

template<typename sph_array_t>
class SphArrayIterator {
public:
  using iterator_category = std::random_access_iterator_tag;
  using value_type = typename sph_array_t::ValueType;
  using difference_type = std::ptrdiff_t;
  using pointer = value_type*;
  using reference = value_type&;
  
  SphArrayIterator(pointer ptr);

  SphArrayIterator& operator++();
  SphArrayIterator operator++(int);

  SphArrayIterator& operator--();
  SphArrayIterator operator--(int);

  SphArrayIterator operator+(const difference_type& index) const;
  SphArrayIterator operator-(const difference_type& index) const;

  difference_type operator-(const SphArrayIterator& it);

  reference operator[](int index);
  pointer operator->();

  reference operator*() const;

  bool operator==(const SphArrayIterator& rhs) const;
  bool operator!=(const SphArrayIterator& rhs) const;
  bool operator<(const SphArrayIterator& rhs) const;
  bool operator<=(const SphArrayIterator& rhs) const;
  bool operator>(const SphArrayIterator& rhs) const;
  bool operator>=(const SphArrayIterator& rhs) const;

  operator SphArrayIterator<typename sph_array_t::ViewTypeConst>() const;

private:
  pointer mPtr;
};

template<typename sph_array_t>
class SphArrayFieldIterator {
public:
  using iterator_category = std::input_iterator_tag;
  using value_type = typename sph_array_t::ValueType;
  using difference_type = std::ptrdiff_t;
  using pointer = value_type*;
  using reference = value_type&;

  using field_type = typename value_type::FieldType;
  using field_pointer = field_type*;
  using field_reference = field_type&;

  SphArrayFieldIterator(pointer ptr);

  field_reference operator*() const;
  field_pointer operator->();

  SphArrayFieldIterator& operator++();
  SphArrayFieldIterator operator++(int);

  bool operator==(const SphArrayFieldIterator& rhs) const;

private:
  pointer mPtr;
};


} //  namespace Spheral

#include "SphArrayInline.hh"

#else

// Forward declare the SphArrayIterator class.
namespace Spheral {
  template<typename sph_array_t> class SphArrayIterator ;
} //  namespace Spheral

#endif //  __Spheral_lvarray_hh__


