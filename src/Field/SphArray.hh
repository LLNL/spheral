#ifndef __Spheral_lvarray_hh__
#define __Spheral_lvarray_hh__

#include "LvArray/Array.hpp"
#include "LvArray/ChaiBuffer.hpp"
#include "chai/ManagedArray.hpp"

namespace Spheral {


template<typename DataType>
class ManagedVector:
  public chai::ManagedArray<DataType>{
  using ManagedArray = chai::ManagedArray<DataType>;

  //friend LvField<DataType>;
public:

  using iterator = DataType*;
  using const_iterator = const DataType*;

  iterator begin() { return ManagedArray::begin(); }
  const_iterator begin() const { return ManagedArray::begin(); }

  iterator end() { return begin() + m_size; }
  const_iterator end() const { return begin() + m_size; }


  RAJA_HOST ManagedVector() : ManagedArray(6, chai::CPU) {}
  RAJA_HOST ManagedVector(size_t elems) : ManagedArray(elems, chai::CPU), m_size(elems) { for (size_t i = 0; i < m_size; i++) ManagedArray::operator[](i) = DataType(); }
  RAJA_HOST ManagedVector(size_t elems, DataType identity) : ManagedArray(elems, chai::CPU), m_size(elems) { for (size_t i = 0; i < m_size; i++) ManagedArray::operator[](i) = identity; }
  //RAJA_HOST ~ManagedVector() { ManagedArray::free(chai::NONE); }
  //RAJA_HOST ~ManagedVector() { destroy(begin(), end()); }

  RAJA_HOST_DEVICE constexpr inline ManagedVector(ManagedVector const& rhs) noexcept : ManagedArray(rhs), m_size(rhs.m_size) {}

  RAJA_HOST void push_back(const DataType& value) {
    if (m_size >= capacity()) ManagedArray::reallocate(capacity() + (capacity() / 2));
    ManagedArray::data()[m_size] = value;
    m_size++;
  }
  RAJA_HOST void push_back(DataType&& value) {
    if (m_size >= capacity()) ManagedArray::reallocate(capacity() + (capacity() / 2));

    ManagedArray::data()[m_size] = std::move(value);
    m_size++;
  }
  template<typename... Args>
  RAJA_HOST
  DataType& emplace_back(Args&&... args) {
    if (m_size >= capacity()) ManagedArray::reallocate(capacity() + (capacity() / 2));

    new(&ManagedArray::data()[m_size]) DataType(std::forward<Args>(args)...);
    return ManagedArray::data()[m_size++];
  }

  RAJA_HOST
  void resize(size_t size) {
    const size_t old_size = m_size;

    if (old_size < size) {
      ManagedArray::reallocate(size);
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

  RAJA_HOST_DEVICE size_t capacity() {return ManagedArray::m_elems;}
  RAJA_HOST_DEVICE size_t size() const {return m_size;}

  RAJA_HOST_DEVICE DataType& operator[](size_t idx) {return ManagedArray::data()[idx]; }
  RAJA_HOST_DEVICE DataType& operator[](size_t idx) const {return ManagedArray::data()[idx]; }

  RAJA_HOST_DEVICE bool operator==(ManagedVector const& rhs) {
    if (m_size != rhs.m_size) return false;
    for (size_t i = 0; i < m_size; i++) if (ManagedArray::data()[i] != rhs[i]) return false;
    return true;
  }


  // *******************************************************
  // Required to Allow ManagedVector to be properly CHAICopyable
  RAJA_HOST_DEVICE ManagedVector<DataType>& operator=(std::nullptr_t) { ManagedArray::operator=(nullptr); return *this; }
  RAJA_HOST_DEVICE void shallowCopy(const ManagedVector& other) {
    m_size=other.m_size;
    ManagedArray::shallowCopy(other);
  }
  // *******************************************************
  
  RAJA_HOST operator std::vector<DataType>() const { return std::vector<DataType>(begin(), end()); }

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


