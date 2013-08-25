//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "peanoHilbertOrderIndicies.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template FieldSpace::FieldList<Dim<1>, KeyTraits::Key> peanoHilbertOrderIndicies<Dim<1> >(const DataBaseSpace::DataBase<Dim<1> >&);
  template FieldSpace::FieldList<Dim<2>, KeyTraits::Key> peanoHilbertOrderIndicies<Dim<2> >(const DataBaseSpace::DataBase<Dim<2> >&);
  template FieldSpace::FieldList<Dim<3>, KeyTraits::Key> peanoHilbertOrderIndicies<Dim<3> >(const DataBaseSpace::DataBase<Dim<3> >&);
}

