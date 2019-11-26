text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "Boundary/Boundary.hh"
#include "CRKSPH/computeCRKSPHSumVolume.cc"

namespace Spheral {
template void computeCRKSPHSumVolume(const ConnectivityMap<Dim< %(ndim)s > >&, 
                                     const TableKernel<Dim< %(ndim)s > >&, 
                                     const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>&,
                                     const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>&,
                                     const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>&,
                                     FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>&);
}
"""
