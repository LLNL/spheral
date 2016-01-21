text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "mortonOrderIndices.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template FieldSpace::FieldList<Dim< %(ndim)s >, KeyTraits::Key> mortonOrderIndices<Dim< %(ndim)s > >(const FieldSpace::FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>&);
  template FieldSpace::FieldList<Dim< %(ndim)s >, KeyTraits::Key> mortonOrderIndices<Dim< %(ndim)s > >(const DataBaseSpace::DataBase<Dim< %(ndim)s > >&);
  template FieldSpace::FieldList<Dim< %(ndim)s >, KeyTraits::Key> mortonOrderIndices<Dim< %(ndim)s > >(const DataBaseSpace::DataBase<Dim< %(ndim)s > >&,
                                                                                     const FieldSpace::FieldList<Dim< %(ndim)s >, int>&);
}

"""
