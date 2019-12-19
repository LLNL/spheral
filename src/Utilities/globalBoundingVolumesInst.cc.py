text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Utilities/globalBoundingVolumes.cc"

namespace Spheral {

  template void globalBoundingBox(const Field<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>& positions,
                                  Dim< %(ndim)s >::Vector& xmin,
                                  Dim< %(ndim)s >::Vector& xmax,
                                  const bool ghost);

  template void globalBoundingBox(const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>& positions,
                                  Dim< %(ndim)s >::Vector& xmin,
                                  Dim< %(ndim)s >::Vector& xmax,
                                  const bool ghost);

  template void globalBoundingVolumes<Dim< %(ndim)s > >(const DataBase<Dim< %(ndim)s > >& dataBase, Dim< %(ndim)s >::ConvexHull& nodeVolume, Dim< %(ndim)s >::ConvexHull& sampleVolume);

}
"""
