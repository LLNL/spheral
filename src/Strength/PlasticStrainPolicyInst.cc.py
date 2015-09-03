text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "PlasticStrainPolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class PlasticStrainPolicy<Dim< %(ndim)s > >;
}

"""
