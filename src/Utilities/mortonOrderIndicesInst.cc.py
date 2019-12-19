text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Utilities/mortonOrderIndices.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template FieldList<Dim< %(ndim)s >, KeyTraits::Key> mortonOrderIndices<Dim< %(ndim)s > >(const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>&);
  template FieldList<Dim< %(ndim)s >, KeyTraits::Key> mortonOrderIndices<Dim< %(ndim)s > >(const DataBase<Dim< %(ndim)s > >&);
  template FieldList<Dim< %(ndim)s >, KeyTraits::Key> mortonOrderIndices<Dim< %(ndim)s > >(const DataBase<Dim< %(ndim)s > >&,
                                                                                           const FieldList<Dim< %(ndim)s >, int>&);
}

"""
