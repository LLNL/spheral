//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Utilities/boundingBox.cc"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template void boundingBox(const vector<Dim<1>::Vector>& positions,
                            Dim<1>::Vector& xmin,
                            Dim<1>::Vector& xmax);

  template void boundingBox(const FieldList<Dim<1>, Dim<1>::Vector>& positions,
                            Dim<1>::Vector& xmin,
                            Dim<1>::Vector& xmax,
                            const bool useGhosts);
#endif

#if defined(SPHERAL_ENABLE_2D)
  template void boundingBox(const vector<Dim<2>::Vector>& positions,
                            Dim<2>::Vector& xmin,
                            Dim<2>::Vector& xmax);

  template void boundingBox(const FieldList<Dim<2>, Dim<2>::Vector>& positions,
                            Dim<2>::Vector& xmin,
                            Dim<2>::Vector& xmax,
                            const bool useGhosts);
#endif

#if defined(SPHERAL_ENABLE_3D)
  template void boundingBox(const vector<Dim<3>::Vector>& positions,
                            Dim<3>::Vector& xmin,
                            Dim<3>::Vector& xmax);

  template void boundingBox(const FieldList<Dim<3>, Dim<3>::Vector>& positions,
                            Dim<3>::Vector& xmin,
                            Dim<3>::Vector& xmax,
                            const bool useGhosts);
#endif
}