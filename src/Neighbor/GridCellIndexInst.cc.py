text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
// Instantiate the static variables.
// Set the maximum index value to 2^21/2, so that (mIndexMax - mIndexMin)^3 will
// be in the range of an unsigned long long (2^64 on the machines I'm using).
#include "Neighbor/GridCellIndex.cc"

namespace Spheral {

  // 2^20
#define INDEXMAX 1048576 - 1
#define INDEXMIN -INDEXMAX - 1
#define XMULTIPLIER (INDEXMAX - INDEXMIN)

  template<> const int GridCellIndex< Dim< %(ndim)s > >::mIndexMax = INDEXMAX; // 2^20 - 1
  template<> const int GridCellIndex< Dim< %(ndim)s > >::mIndexMin = INDEXMIN;

}
"""
