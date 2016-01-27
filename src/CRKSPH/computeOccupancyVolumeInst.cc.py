text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "computeOccupancyVolume.cc"

namespace Spheral {
  namespace CRKSPHSpace {
    template void computeOccupancyVolume(const ConnectivityMap<Dim< %(ndim)s > >&, 
                                         const TableKernel<Dim< %(ndim)s > >&, 
                                         const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>&,
                                         const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>&,
                                         FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>&);
  }
}
"""
