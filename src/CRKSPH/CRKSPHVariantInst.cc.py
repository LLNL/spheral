text = """
//------------------------------------------------------------------------------
// Explict instantiation.
//------------------------------------------------------------------------------
#include "CRKSPH/CRKSPHVariant.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  namespace CRKSPHSpace {
    template class CRKSPHVariant< Dim< %(ndim)s > >;
  }
}
"""
