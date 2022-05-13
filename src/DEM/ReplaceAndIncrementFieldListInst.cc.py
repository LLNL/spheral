text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "DEM/ReplaceAndIncrementFieldList.cc"

namespace Spheral {
  template class ReplaceAndIncrementFieldList<Dim< %(ndim)s >,std::vector<Dim< %(ndim)s>::Vector>>;
}
"""
