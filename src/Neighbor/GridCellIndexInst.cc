//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Neighbor/GridCellIndex.cc"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)

  // 2^20
#define INDEXMAX 1048576 - 1
#define INDEXMIN -INDEXMAX - 1
#define XMULTIPLIER (INDEXMAX - INDEXMIN)

  template<> const int GridCellIndex< Dim<1> >::mIndexMax = INDEXMAX; // 2^20 - 1
  template<> const int GridCellIndex< Dim<1> >::mIndexMin = INDEXMIN;

#endif

#if defined(SPHERAL_ENABLE_2D)

  // 2^20
#define INDEXMAX 1048576 - 1
#define INDEXMIN -INDEXMAX - 1
#define XMULTIPLIER (INDEXMAX - INDEXMIN)

  template<> const int GridCellIndex< Dim<2> >::mIndexMax = INDEXMAX; // 2^20 - 1
  template<> const int GridCellIndex< Dim<2> >::mIndexMin = INDEXMIN;

#endif

#if defined(SPHERAL_ENABLE_3D)

  // 2^20
#define INDEXMAX 1048576 - 1
#define INDEXMIN -INDEXMAX - 1
#define XMULTIPLIER (INDEXMAX - INDEXMIN)

  template<> const int GridCellIndex< Dim<3> >::mIndexMax = INDEXMAX; // 2^20 - 1
  template<> const int GridCellIndex< Dim<3> >::mIndexMin = INDEXMIN;

#endif
}