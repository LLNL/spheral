//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Utilities/globalBoundingVolumes.cc"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)

  template void globalBoundingBox(const Field<Dim<1>, Dim<1>::Vector>& positions,
                                  Dim<1>::Vector& xmin,
                                  Dim<1>::Vector& xmax,
                                  const bool ghost);

  template void globalBoundingBox(const FieldList<Dim<1>, Dim<1>::Vector>& positions,
                                  Dim<1>::Vector& xmin,
                                  Dim<1>::Vector& xmax,
                                  const bool ghost);

  template void globalBoundingVolumes<Dim<1> >(const DataBase<Dim<1> >& dataBase, Dim<1>::ConvexHull& nodeVolume, Dim<1>::ConvexHull& sampleVolume);

#endif

#if defined(SPHERAL_ENABLE_2D)

  template void globalBoundingBox(const Field<Dim<2>, Dim<2>::Vector>& positions,
                                  Dim<2>::Vector& xmin,
                                  Dim<2>::Vector& xmax,
                                  const bool ghost);

  template void globalBoundingBox(const FieldList<Dim<2>, Dim<2>::Vector>& positions,
                                  Dim<2>::Vector& xmin,
                                  Dim<2>::Vector& xmax,
                                  const bool ghost);

  template void globalBoundingVolumes<Dim<2> >(const DataBase<Dim<2> >& dataBase, Dim<2>::ConvexHull& nodeVolume, Dim<2>::ConvexHull& sampleVolume);

#endif

#if defined(SPHERAL_ENABLE_3D)

  template void globalBoundingBox(const Field<Dim<3>, Dim<3>::Vector>& positions,
                                  Dim<3>::Vector& xmin,
                                  Dim<3>::Vector& xmax,
                                  const bool ghost);

  template void globalBoundingBox(const FieldList<Dim<3>, Dim<3>::Vector>& positions,
                                  Dim<3>::Vector& xmin,
                                  Dim<3>::Vector& xmax,
                                  const bool ghost);

  template void globalBoundingVolumes<Dim<3> >(const DataBase<Dim<3> >& dataBase, Dim<3>::ConvexHull& nodeVolume, Dim<3>::ConvexHull& sampleVolume);

#endif
}