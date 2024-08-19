//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Utilities/mortonOrderIndices.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template FieldList<Dim<1>, KeyTraits::Key> mortonOrderIndices<Dim<1> >(const FieldList<Dim<1>, Dim<1>::Vector>&);
  template FieldList<Dim<1>, KeyTraits::Key> mortonOrderIndices<Dim<1> >(const DataBase<Dim<1> >&);
  template FieldList<Dim<1>, KeyTraits::Key> mortonOrderIndices<Dim<1> >(const DataBase<Dim<1> >&,
                                                                                           const FieldList<Dim<1>, int>&);
#endif

#if defined(SPHERAL_ENABLE_2D)
  template FieldList<Dim<2>, KeyTraits::Key> mortonOrderIndices<Dim<2> >(const FieldList<Dim<2>, Dim<2>::Vector>&);
  template FieldList<Dim<2>, KeyTraits::Key> mortonOrderIndices<Dim<2> >(const DataBase<Dim<2> >&);
  template FieldList<Dim<2>, KeyTraits::Key> mortonOrderIndices<Dim<2> >(const DataBase<Dim<2> >&,
                                                                                           const FieldList<Dim<2>, int>&);
#endif

#if defined(SPHERAL_ENABLE_3D)
  template FieldList<Dim<3>, KeyTraits::Key> mortonOrderIndices<Dim<3> >(const FieldList<Dim<3>, Dim<3>::Vector>&);
  template FieldList<Dim<3>, KeyTraits::Key> mortonOrderIndices<Dim<3> >(const DataBase<Dim<3> >&);
  template FieldList<Dim<3>, KeyTraits::Key> mortonOrderIndices<Dim<3> >(const DataBase<Dim<3> >&,
                                                                                           const FieldList<Dim<3>, int>&);
#endif
}