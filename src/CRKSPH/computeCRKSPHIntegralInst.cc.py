text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "CRKSPH/computeCRKSPHIntegral.cc"

namespace Spheral {
template std::pair<Dim< %(ndim)s >::Vector,Dim< %(ndim)s >::Vector> computeCRKSPHIntegral(const ConnectivityMap<Dim< %(ndim)s > >&, 
                                                                                          const TableKernel<Dim< %(ndim)s > >&, 
                                                                                          const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>&,
                                                                                          const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>&,
                                                                                          const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>&,
                                                                                          size_t, const int, size_t, const int, int, const int,
                                                                                          Dim< %(ndim)s >::Vector, Dim< %(ndim)s >::Vector);
}

"""
