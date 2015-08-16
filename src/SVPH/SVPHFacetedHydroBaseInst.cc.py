text = """
//------------------------------------------------------------------------------
// Explict instantiation.
//------------------------------------------------------------------------------
#include "SVPHFacetedHydroBase.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  namespace SVPHSpace {
    template class SVPHFacetedHydroBase< Dim< %(ndim)s > >;
  }
}
"""
