//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "BulkModulusPolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class BulkModulusPolicy<Dim<1> >;
  template class BulkModulusPolicy<Dim<2> >;
  template class BulkModulusPolicy<Dim<3> >;
}
