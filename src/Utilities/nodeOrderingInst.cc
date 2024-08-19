//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Utilities/nodeOrdering.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template FieldList<Dim<1>, int> nodeOrdering<Dim<1> >(const FieldList<Dim<1>, KeyTraits::Key>&);
#endif

#if defined(SPHERAL_ENABLE_2D)
  template FieldList<Dim<2>, int> nodeOrdering<Dim<2> >(const FieldList<Dim<2>, KeyTraits::Key>&);
#endif

#if defined(SPHERAL_ENABLE_3D)
  template FieldList<Dim<3>, int> nodeOrdering<Dim<3> >(const FieldList<Dim<3>, KeyTraits::Key>&);
#endif
}