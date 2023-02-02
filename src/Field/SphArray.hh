#ifndef __Spheral_lvarray_hh__
#define __Spheral_lvarray_hh__

#include "LvArray/Array.hpp"
#include "LvArray/ChaiBuffer.hpp"

namespace Spheral {

template<typename DATA_TYPE>
using SphArray = LvArray::Array<DATA_TYPE, 1, camp::idx_seq<0>, std::ptrdiff_t, LvArray::ChaiBuffer>;

template<typename DATA_TYPE>
using SphArrayView = LvArray::ArrayView<DATA_TYPE, 1, 0, std::ptrdiff_t, LvArray::ChaiBuffer>;

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


};

#include "SphArrayInline.hh"

#else

// Forward declare the SphArrayIterator class.
namespace Spheral {
  template<typename sph_array_t> class SphArrayIterator ;
} //  namespace Spheral

#endif //  __Spheral_lvarray_hh__


