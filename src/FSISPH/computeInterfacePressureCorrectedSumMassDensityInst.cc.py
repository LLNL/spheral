text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "FSISPH/computeInterfacePressureCorrectedSumMassDensity.cc"
#include "Geometry/Dimension.hh"
namespace Spheral {
  template void computeInterfacePressureCorrectedSumMassDensity(const ConnectivityMap<Dim< %(ndim)s > >&, 
                                                                const TableKernel<Dim< %(ndim)s > >&, 
                                                                const std::vector<int>&,
                                                                const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>&,
                                                                const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>&,
                                                                const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>&,
                                                                const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>&,
                                                                const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>&,
                                                                const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>&,
                                                                      FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>&);
}
"""