//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "mortonOrderIndices.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template FieldSpace::FieldList<Dim<1>, KeyTraits::Key> mortonOrderIndices<Dim<1> >(const FieldSpace::FieldList<Dim<1>, Dim<1>::Vector>&);
  template FieldSpace::FieldList<Dim<2>, KeyTraits::Key> mortonOrderIndices<Dim<2> >(const FieldSpace::FieldList<Dim<2>, Dim<2>::Vector>&);
  template FieldSpace::FieldList<Dim<3>, KeyTraits::Key> mortonOrderIndices<Dim<3> >(const FieldSpace::FieldList<Dim<3>, Dim<3>::Vector>&);

  template FieldSpace::FieldList<Dim<1>, KeyTraits::Key> mortonOrderIndices<Dim<1> >(const DataBaseSpace::DataBase<Dim<1> >&);
  template FieldSpace::FieldList<Dim<2>, KeyTraits::Key> mortonOrderIndices<Dim<2> >(const DataBaseSpace::DataBase<Dim<2> >&);
  template FieldSpace::FieldList<Dim<3>, KeyTraits::Key> mortonOrderIndices<Dim<3> >(const DataBaseSpace::DataBase<Dim<3> >&);

  template FieldSpace::FieldList<Dim<1>, KeyTraits::Key> mortonOrderIndices<Dim<1> >(const DataBaseSpace::DataBase<Dim<1> >&,
                                                                                     const FieldSpace::FieldList<Dim<1>, int>&);
  template FieldSpace::FieldList<Dim<2>, KeyTraits::Key> mortonOrderIndices<Dim<2> >(const DataBaseSpace::DataBase<Dim<2> >&,
                                                                                     const FieldSpace::FieldList<Dim<2>, int>&);
  template FieldSpace::FieldList<Dim<3>, KeyTraits::Key> mortonOrderIndices<Dim<3> >(const DataBaseSpace::DataBase<Dim<3> >&,
                                                                                     const FieldSpace::FieldList<Dim<3>, int>&);
}

