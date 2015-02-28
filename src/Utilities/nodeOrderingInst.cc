//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "nodeOrdering.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template FieldSpace::FieldList<Dim<1>, int> nodeOrdering<Dim<1> >(const FieldSpace::FieldList<Dim<1>, KeyTraits::Key>&);
  template FieldSpace::FieldList<Dim<2>, int> nodeOrdering<Dim<2> >(const FieldSpace::FieldList<Dim<2>, KeyTraits::Key>&);
  template FieldSpace::FieldList<Dim<3>, int> nodeOrdering<Dim<3> >(const FieldSpace::FieldList<Dim<3>, KeyTraits::Key>&);
}
