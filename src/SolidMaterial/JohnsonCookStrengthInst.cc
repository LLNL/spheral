//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "JohnsonCookStrength.cc"
#include "Geometry/Dimension.hh"

using namespace Spheral::Material;

namespace Spheral {
  namespace SolidMaterial {
    template class JohnsonCookStrength<Dim<1> >;
    template class JohnsonCookStrength<Dim<2> >;
    template class JohnsonCookStrength<Dim<3> >;
  }
}
