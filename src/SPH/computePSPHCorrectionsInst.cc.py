text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "computePSPHCorrections.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  namespace SPHSpace {
    template void computePSPHCorrections(const ConnectivityMap<Dim< %(ndim)s > >&, 
                                                 const TableKernel<Dim< %(ndim)s > >&, 
                                                 const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>&,
                                                 const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>&,
                                                 const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>&,
                                                 const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>&,
                                                 const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>&,
                                                 const bool,
                                                 FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>&,
                                                 FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>&,
                                                 FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>&,
                                                 FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>&);
  }
}

"""
