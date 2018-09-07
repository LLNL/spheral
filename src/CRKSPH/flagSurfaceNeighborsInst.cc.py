text = """
//------------------------------------------------------------------------------
// Explict instantiation.
//------------------------------------------------------------------------------
#include "flagSurfaceNeighbors.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
template void flagSurfaceNeighbors<Dim<%(ndim)s>>(FieldList<Dim<%(ndim)s>, int>& surfacePoint,
                                                  const ConnectivityMap<Dim<%(ndim)s>>& connectivityMap);
}
"""
