text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "SPH/computeSPHSumMassDensity.cc"
#include "Geometry/Dimension.hh"
#include "Kernel/TableKernel.hh"

namespace Spheral {
  template void computeSPHSumMassDensity(const ConnectivityMap<%(Dimension)s>&, 
                                         const TableKernel<%(Dimension)s>&, 
                                         const bool,
                                         const FieldList<%(Dimension)s, %(Dimension)s::Vector>&,
                                         const FieldList<Dim<%(ndim)s, %(Dimension)s::Scalar>&,
                                         const FieldList<Dim<%(ndim)s >, %(Dimension)s::SymTensor>&,
                                         FieldList<%(Dimension)s, %(Dimension)s::Scalar>&);
}
"""
