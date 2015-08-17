text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "computeSPHOmegaGradhCorrection.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  namespace SPHSpace {
    template void computeSPHOmegaGradhCorrection(const ConnectivityMap<Dim< %(ndim)s > >&, 
                                                 const TableKernel<Dim< %(ndim)s > >&, 
                                                 const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>&,
                                                 const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>&,
                                                 FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>&);
  }
}

"""
