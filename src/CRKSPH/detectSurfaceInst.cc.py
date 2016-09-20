text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "computeCRKSPHMoments.cc"

namespace Spheral {
namespace CRKSPHSpace {

    template void detectSurface(const ConnectivityMap< Dim< %(ndim)s > >& connectivityMap,
                                const FieldList< Dim< %(ndim)s >,  Dim< %(ndim)s >::Scalar>& m0,
                                const FieldList< Dim< %(ndim)s >,  Dim< %(ndim)s >::Vector>& m1,
                                const double detectThreshold,
                                const double detectRange,
                                const double sweepAngle,
                                FieldList< Dim< %(ndim)s >,  Dim< %(ndim)s >::Scalar>& surfNorm);
}
}
"""