//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "copy2polytope.cc"

namespace Spheral {
#if defined(SPHERAL_ENABLE_2D)
  template void copy2polytope(const FieldList<Dim<2>, Dim<2>::FacetedVolume>& cells,
                              polytope::Tessellation<2, double>& mesh);
#endif

#if defined(SPHERAL_ENABLE_3D)
  template void copy2polytope(const FieldList<Dim<3>, Dim<3>::FacetedVolume>& cells,
                              polytope::Tessellation<3, double>& mesh);
#endif
}
