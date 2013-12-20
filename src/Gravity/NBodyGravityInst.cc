//------------------------------------------------------------------------------
// Explict instantiations.
//------------------------------------------------------------------------------
#include "NBodyGravity.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  namespace GravitySpace {

    template class NBodyGravity<Dim<1> >;
    template class NBodyGravity<Dim<2> >;
    template class NBodyGravity<Dim<3> >;

  } // end namespace GravitySpace
} // end namespace Spheral

