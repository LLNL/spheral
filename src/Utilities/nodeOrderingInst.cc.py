text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "nodeOrdering.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template FieldSpace::FieldList<Dim< %(ndim)s >, int> nodeOrdering<Dim< %(ndim)s > >(const FieldSpace::FieldList<Dim< %(ndim)s >, KeyTraits::Key>&);
}
"""
