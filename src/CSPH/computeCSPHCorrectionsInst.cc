//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "computeCSPHCorrections.cc"

namespace Spheral {
  namespace CSPHSpace {
    template void computeCSPHCorrections(const ConnectivityMap<Dim<1> >&, 
                                         const TableKernel<Dim<1> >&, 
                                         const FieldList<Dim<1>, Dim<1>::Scalar>&,
                                         const FieldList<Dim<1>, Dim<1>::Vector>&,
                                         const FieldList<Dim<1>, Dim<1>::SymTensor>&,
                                         const bool,
                                         FieldList<Dim<1>, Dim<1>::Scalar>&,
                                         FieldList<Dim<1>, Dim<1>::Vector>&,
                                         FieldList<Dim<1>, Dim<1>::SymTensor>&,
                                         FieldList<Dim<1>, Dim<1>::Scalar>&,
                                         FieldList<Dim<1>, Dim<1>::Scalar>&,
                                         FieldList<Dim<1>, Dim<1>::Vector>&,
                                         FieldList<Dim<1>, Dim<1>::Vector>&,
                                         FieldList<Dim<1>, Dim<1>::Tensor>&,
                                         FieldList<Dim<1>, Dim<1>::Vector>&,
                                         FieldList<Dim<1>, Dim<1>::Tensor>&);
    template void computeCSPHCorrections(const ConnectivityMap<Dim<2> >&, 
                                         const TableKernel<Dim<2> >&, 
                                         const FieldList<Dim<2>, Dim<2>::Scalar>&,
                                         const FieldList<Dim<2>, Dim<2>::Vector>&,
                                         const FieldList<Dim<2>, Dim<2>::SymTensor>&,
                                         const bool,
                                         FieldList<Dim<2>, Dim<2>::Scalar>&,
                                         FieldList<Dim<2>, Dim<2>::Vector>&,
                                         FieldList<Dim<2>, Dim<2>::SymTensor>&,
                                         FieldList<Dim<2>, Dim<2>::Scalar>&,
                                         FieldList<Dim<2>, Dim<2>::Scalar>&,
                                         FieldList<Dim<2>, Dim<2>::Vector>&,
                                         FieldList<Dim<2>, Dim<2>::Vector>&,
                                         FieldList<Dim<2>, Dim<2>::Tensor>&,
                                         FieldList<Dim<2>, Dim<2>::Vector>&,
                                         FieldList<Dim<2>, Dim<2>::Tensor>&);
    template void computeCSPHCorrections(const ConnectivityMap<Dim<3> >&, 
                                         const TableKernel<Dim<3> >&, 
                                         const FieldList<Dim<3>, Dim<3>::Scalar>&,
                                         const FieldList<Dim<3>, Dim<3>::Vector>&,
                                         const FieldList<Dim<3>, Dim<3>::SymTensor>&,
                                         const bool,
                                         FieldList<Dim<3>, Dim<3>::Scalar>&,
                                         FieldList<Dim<3>, Dim<3>::Vector>&,
                                         FieldList<Dim<3>, Dim<3>::SymTensor>&,
                                         FieldList<Dim<3>, Dim<3>::Scalar>&,
                                         FieldList<Dim<3>, Dim<3>::Scalar>&,
                                         FieldList<Dim<3>, Dim<3>::Vector>&,
                                         FieldList<Dim<3>, Dim<3>::Vector>&,
                                         FieldList<Dim<3>, Dim<3>::Tensor>&,
                                         FieldList<Dim<3>, Dim<3>::Vector>&,
                                         FieldList<Dim<3>, Dim<3>::Tensor>&);
  }
}

