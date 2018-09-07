text = """
//------------------------------------------------------------------------------
// Explict instantiation.
//------------------------------------------------------------------------------
#include "SVPHFacetedHydroBase.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class SVPHFacetedHydroBase< Dim< %(ndim)s > >;
}
"""
