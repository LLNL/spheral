text = """
//------------------------------------------------------------------------------
// Explict instantiation.
//------------------------------------------------------------------------------
#include "SVPH/SVPHFacetedHydroBase.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  namespace SVPHSpace {
    template class SVPHFacetedHydroBase< Dim< %(ndim)s > >;
  }
}
"""
