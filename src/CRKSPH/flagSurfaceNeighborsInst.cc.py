text = """
//------------------------------------------------------------------------------
// Explict instantiation.
//------------------------------------------------------------------------------
#include "CRKSPH/flagSurfaceNeighbors.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  namespace CRKSPHSpace {
    template void flagSurfaceNeighbors<Dim<%(ndim)s>>(FieldSpace::FieldList<Dim<%(ndim)s>, int>& surfacePoint,
                                                      const NeighborSpace::ConnectivityMap<Dim<%(ndim)s>>& connectivityMap);
  }
}
"""
