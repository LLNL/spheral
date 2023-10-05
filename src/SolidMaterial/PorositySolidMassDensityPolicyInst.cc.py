text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "SolidMaterial/PorositySolidMassDensityPolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class PorositySolidMassDensityPolicy<Dim< %(ndim)s > >;
}

"""
