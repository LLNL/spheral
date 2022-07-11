text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "SPH/computeSPHSumMassDensity.cc"
#include "Geometry/Dimension.hh"
#include "Kernel/TableKernel.hh"

namespace Spheral {
  template void computeSPHSumMassDensity(const ConnectivityMap<%(Dim)s>&, 
                                         const TableKernel<%(Dim)s>&, 
                                         const bool,
                                         const FieldList<%(Dim)s, %(Dim)s::Vector>&,
                                         const FieldList<%(Dim)s, %(Dim)s::Scalar>&,
                                         const FieldList<%(Dim)s, %(Dim)s::SymTensor>&,
                                         FieldList<%(Dim)s, %(Dim)s::Scalar>&);
}
"""
