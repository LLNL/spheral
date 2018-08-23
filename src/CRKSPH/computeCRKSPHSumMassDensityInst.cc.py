text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "Boundary/Boundary.hh"

#ifdef _OPENMP
#include "CRKSPH/computeCRKSPHSumMassDensity_OpenMP.cc"
#else
#include "CRKSPH/computeCRKSPHSumMassDensity.cc"
#endif

namespace Spheral {
  namespace CRKSPHSpace {
    template void computeCRKSPHSumMassDensity(const ConnectivityMap<Dim< %(ndim)s > >&, 
                                              const TableKernel<Dim< %(ndim)s > >&, 
                                              const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>& position,
                                              const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& mass,
                                              const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& vol,
                                              const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>& H,
                                              const FieldList<Dim< %(ndim)s >, int>& voidPoint,
                                              FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& massDensity);
  }
}
"""
