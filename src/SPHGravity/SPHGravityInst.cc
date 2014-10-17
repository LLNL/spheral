//------------------------------------------------------------------------------
// Explict instantiations.
//------------------------------------------------------------------------------
#include "SPHGravity.cc"

namespace Spheral {
  namespace GravitySpace {

    template class SPHGravity<Dim<1> >;
    template class SPHGravity<Dim<2> >;
    template class SPHGravity<Dim<3> >;

  } // end namespace GravitySpace
} // end namespace Spheral

