//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "computeCRKSPHIntegral.cc"

namespace Spheral {
  namespace CRKSPHSpace {
    template std::pair<Dim<1>::Vector,Dim<1>::Vector> computeCRKSPHIntegral(const ConnectivityMap<Dim<1> >&, 
                                                                          const TableKernel<Dim<1> >&, 
                                                                          const FieldList<Dim<1>, Dim<1>::Scalar>&,
                                                                          const FieldList<Dim<1>, Dim<1>::Vector>&,
                                                                          const FieldList<Dim<1>, Dim<1>::SymTensor>&,
                                                                          size_t, const int, size_t, const int, int, const int,
                                                                          Dim<1>::Vector, Dim<1>::Vector);
    template std::pair<Dim<2>::Vector,Dim<2>::Vector> computeCRKSPHIntegral(const ConnectivityMap<Dim<2> >&, 
                                                                          const TableKernel<Dim<2> >&, 
                                                                          const FieldList<Dim<2>, Dim<2>::Scalar>&,
                                                                          const FieldList<Dim<2>, Dim<2>::Vector>&,
                                                                          const FieldList<Dim<2>, Dim<2>::SymTensor>&,
                                                                          size_t, const int, size_t, const int, int, const int,
                                                                          Dim<2>::Vector, Dim<2>::Vector);
    template std::pair<Dim<3>::Vector,Dim<3>::Vector> computeCRKSPHIntegral(const ConnectivityMap<Dim<3> >&, 
                                                                          const TableKernel<Dim<3> >&, 
                                                                          const FieldList<Dim<3>, Dim<3>::Scalar>&,
                                                                          const FieldList<Dim<3>, Dim<3>::Vector>&,
                                                                          const FieldList<Dim<3>, Dim<3>::SymTensor>&,
                                                                          size_t, const int, size_t, const int, int, const int,
                                                                          Dim<3>::Vector, Dim<3>::Vector);
  }
}

