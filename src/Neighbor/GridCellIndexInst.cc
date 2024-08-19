//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Neighbor/GridCellIndex.cc"

#define INDEXMAX 1048576 - 1 // 2^20 - 1
#define INDEXMIN -INDEXMAX - 1 // -2^20

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template<> const int GridCellIndex< Dim<1> >::mIndexMax = INDEXMAX;
  template<> const int GridCellIndex< Dim<1> >::mIndexMin = INDEXMIN;
#endif

#if defined(SPHERAL_ENABLE_2D)
  template<> const int GridCellIndex< Dim<2> >::mIndexMax = INDEXMAX;
  template<> const int GridCellIndex< Dim<2> >::mIndexMin = INDEXMIN;
#endif

#if defined(SPHERAL_ENABLE_3D)
  template<> const int GridCellIndex< Dim<3> >::mIndexMax = INDEXMAX;
  template<> const int GridCellIndex< Dim<3> >::mIndexMin = INDEXMIN;
#endif
}
