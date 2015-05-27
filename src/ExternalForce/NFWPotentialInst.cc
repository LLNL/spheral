//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "NFWPotential.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  namespace PhysicsSpace {
    using namespace Material;
    template class NFWPotential<Dim<1> >;
    template class NFWPotential<Dim<2> >;
    template class NFWPotential<Dim<3> >;
  }
}
