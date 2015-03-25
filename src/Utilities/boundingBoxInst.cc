//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "boundingBox.cc"

namespace Spheral {
  template void boundingBox(const vector<Dim<1>::Vector>& positions,
                            Dim<1>::Vector& xmin,
                            Dim<1>::Vector& xmax);
  template void boundingBox(const vector<Dim<2>::Vector>& positions,
                            Dim<2>::Vector& xmin,
                            Dim<2>::Vector& xmax);
  template void boundingBox(const vector<Dim<3>::Vector>& positions,
                            Dim<3>::Vector& xmin,
                            Dim<3>::Vector& xmax);

  template void boundingBox(const FieldSpace::FieldList<Dim<1>, Dim<1>::Vector>& positions,
                            Dim<1>::Vector& xmin,
                            Dim<1>::Vector& xmax,
                            const bool useGhosts);
  template void boundingBox(const FieldSpace::FieldList<Dim<2>, Dim<2>::Vector>& positions,
                            Dim<2>::Vector& xmin,
                            Dim<2>::Vector& xmax,
                            const bool useGhosts);
  template void boundingBox(const FieldSpace::FieldList<Dim<3>, Dim<3>::Vector>& positions,
                            Dim<3>::Vector& xmin,
                            Dim<3>::Vector& xmax,
                            const bool useGhosts);
}
