text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "CRKSPH/detectSurface.cc"

namespace Spheral {
template void detectSurface(const ConnectivityMap< Dim< %(ndim)s > >& connectivityMap,
                            const FieldList< Dim< %(ndim)s >,  Dim< %(ndim)s >::Scalar>& m0,
                            const FieldList< Dim< %(ndim)s >,  Dim< %(ndim)s >::Vector>& m1,
                            const FieldList< Dim< %(ndim)s >,  Dim< %(ndim)s >::Vector>& position,
                            const FieldList< Dim< %(ndim)s >,  Dim< %(ndim)s >::SymTensor>& H,
                            const double detectThreshold,
                            const double detectRange,
                            const double sweepAngle,
                            FieldList< Dim< %(ndim)s >, int>& surfacePoint);
}
"""
