text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "computeCRKSPHEvaluation.cc"

namespace Spheral {
  namespace CRKSPHSpace {
    template void computeCRKSPHEvaluation(const ConnectivityMap<Dim< %(ndim)s > >&, 
                                         const TableKernel<Dim< %(ndim)s > >&, 
                                         const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>&,
                                         const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>&,
                                         const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>&,
                                         size_t, const int, Dim< %(ndim)s >::Vector,
                                         const bool, Dim< %(ndim)s >::Scalar&, Dim< %(ndim)s >::Vector&);
  }
}

"""
