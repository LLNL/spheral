text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Strength/PlasticStrainPolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class PlasticStrainPolicy<Dim< %(ndim)s > >;
}

"""
