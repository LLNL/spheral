#ifndef __Spheral_lvarray_hh__
#define __Spheral_lvarray_hh__

#include "LvArray/Array.hpp"
#include "LvArray/ChaiBuffer.hpp"

namespace Spheral {

template<typename DATA_TYPE>
using SphArray = LvArray::Array<DATA_TYPE, 1, camp::idx_seq<0>, std::ptrdiff_t, LvArray::ChaiBuffer>;

template<typename DATA_TYPE>
using SphArrayView = LvArray::ArrayView<DATA_TYPE, 1, 0, std::ptrdiff_t, LvArray::ChaiBuffer>;

}

#endif //  __Spheral_lvarray_hh__


