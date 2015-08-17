text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "Boundary/Boundary.hh"
#include "SolidSPH/NodeCoupling.hh"
#include "computeSolidCRKSPHSumMassDensity.cc"

namespace Spheral {
  namespace CRKSPHSpace {
    template void computeSolidCRKSPHSumMassDensity(const ConnectivityMap<Dim< %(ndim)s > >&, 
                                                   const TableKernel<Dim< %(ndim)s > >&, 
                                                   const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>&,
                                                   const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>&,
                                                   const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>&,
                                                   const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>&,
                                                   const NodeCoupling&,
                                                   FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>&);
  }
}
"""
