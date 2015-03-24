//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "SteinbergGuinanStrength.cc"
#include "Geometry/Dimension.hh"

using namespace Spheral::Material;

namespace Spheral {
  namespace SolidMaterial {
    template class SteinbergGuinanStrength<Dim<1> >;
    template class SteinbergGuinanStrength<Dim<2> >;
    template class SteinbergGuinanStrength<Dim<3> >;
  }
}
