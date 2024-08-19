//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "DataBase/DataBase.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template class DataBase< Dim<1> >;
#endif

#if defined(SPHERAL_ENABLE_2D)
  template class DataBase< Dim<2> >;
#endif

#if defined(SPHERAL_ENABLE_3D)
  template class DataBase< Dim<3> >;
#endif
}