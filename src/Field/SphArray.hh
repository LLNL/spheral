#ifndef __Spheral_lvarray_hh__
#define __Spheral_lvarray_hh__

#include "LvArray/Array.hpp"
#include "LvArray/ChaiBuffer.hpp"

namespace Spheral {

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


