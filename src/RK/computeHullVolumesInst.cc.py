text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "RK/computeHullVolumes.cc"

namespace Spheral {
  template void computeHullVolumes(const ConnectivityMap<Dim< %(ndim)s > >&, 
                                   const Dim< %(ndim)s >::Scalar kernelExtent,
                                   const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>&,
                                   const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>&,
                                   FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>&);
}

"""
