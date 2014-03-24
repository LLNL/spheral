//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "peanoHilbertOrderIndices.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template FieldSpace::FieldList<Dim<1>, KeyTraits::Key> peanoHilbertOrderIndices<Dim<1> >(const DataBaseSpace::DataBase<Dim<1> >&);
  template FieldSpace::FieldList<Dim<2>, KeyTraits::Key> peanoHilbertOrderIndices<Dim<2> >(const DataBaseSpace::DataBase<Dim<2> >&);
  template FieldSpace::FieldList<Dim<3>, KeyTraits::Key> peanoHilbertOrderIndices<Dim<3> >(const DataBaseSpace::DataBase<Dim<3> >&);
}

