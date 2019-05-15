text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "editMultimaterialSurfaceTopology.cc"

namespace Spheral {
template void editMultimaterialSurfaceTopology(FieldList< Dim< %(ndim)s >, int>& surfacePoint,
                                               ConnectivityMap< Dim< %(ndim)s > >& connectivityMap);
}
"""
