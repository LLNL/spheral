text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "detectSurface.cc"

namespace Spheral {
namespace CRKSPHSpace {

    template void detectSurface(const NeighborSpace::ConnectivityMap< Dim< %(ndim)s > >& connectivityMap,
                                const FieldSpace::FieldList< Dim< %(ndim)s >,  Dim< %(ndim)s >::Scalar>& m0,
                                const FieldSpace::FieldList< Dim< %(ndim)s >,  Dim< %(ndim)s >::Vector>& m1,
                                const FieldSpace::FieldList< Dim< %(ndim)s >,  Dim< %(ndim)s >::Vector>& position,
                                const FieldSpace::FieldList< Dim< %(ndim)s >,  Dim< %(ndim)s >::SymTensor>& H,
                                const double detectThreshold,
                                const double detectRange,
                                const double sweepAngle,
                                FieldSpace::FieldList< Dim< %(ndim)s >, int>& surfacePoint);
}
}
"""
