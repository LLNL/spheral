text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "Boundary/Boundary.hh"
#include "CRKSPH/computeHullSumMassDensity.cc"
#include "SolidSPH/DamagedNodeCouplingWithFrags.hh"

namespace Spheral {
template void computeHullSumMassDensity(const ConnectivityMap<Dim< %(ndim)s > >&, 
                                        const TableKernel<Dim< %(ndim)s > >&,
                                        const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>&,
                                        const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>&,
                                        const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>&,
                                        const NodeCoupling&,
                                        FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>&);
}
"""
