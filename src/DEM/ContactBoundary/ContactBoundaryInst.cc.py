text = """
//------------------------------------------------------------------------------
// Explict instantiation.
//------------------------------------------------------------------------------
#include "DEM/ContactBoundary/ContactBoundary.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class ContactBoundary< Dim< %(ndim)s > >;
}
"""
