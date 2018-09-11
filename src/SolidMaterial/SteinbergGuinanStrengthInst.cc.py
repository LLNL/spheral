text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "SolidMaterial/SteinbergGuinanStrength.cc"
#include "Geometry/Dimension.hh"

using namespace Spheral::Material;

namespace Spheral {
  namespace SolidMaterial {
    template class SteinbergGuinanStrength<Dim< %(ndim)s > >;
  }
}
"""
