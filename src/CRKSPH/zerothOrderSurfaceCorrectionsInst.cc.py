text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "CRKSPH/zerothOrderSurfaceCorrections.cc"

namespace Spheral {
template void
zerothOrderSurfaceCorrections(FieldList<Dim<%(ndim)s>, Dim<%(ndim)s>::Scalar>& A,
                              FieldList<Dim<%(ndim)s>, Dim<%(ndim)s>::Vector>& B,
                              FieldList<Dim<%(ndim)s>, Dim<%(ndim)s>::Tensor>& C,
                              FieldList<Dim<%(ndim)s>, Dim<%(ndim)s>::Vector>& gradA,
                              FieldList<Dim<%(ndim)s>, Dim<%(ndim)s>::Tensor>& gradB,
                              FieldList<Dim<%(ndim)s>, Dim<%(ndim)s>::ThirdRankTensor>& gradC,
                              const FieldList<Dim<%(ndim)s>, Dim<%(ndim)s>::Scalar>& m0,
                              const FieldList<Dim<%(ndim)s>, Dim<%(ndim)s>::Vector>& gradm0,
                              const FieldList<Dim<%(ndim)s>, int>& surfacePoint);
}
"""
