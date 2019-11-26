text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "Boundary/Boundary.hh"
#include "CRKSPH/computeCRKSPHSumMassDensity.cc"

namespace Spheral {
template void computeCRKSPHSumMassDensity(const ConnectivityMap<Dim< %(ndim)s > >&, 
                                          const TableKernel<Dim< %(ndim)s > >&, 
                                          const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>& position,
                                          const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& mass,
                                          const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& vol,
                                          const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>& H,
                                          FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& massDensity);
}
"""
