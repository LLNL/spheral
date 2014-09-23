//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "Boundary/Boundary.hh"
#include "computeCSPHSumMassDensity.cc"

namespace Spheral {
  namespace CSPHSpace {
    template void computeCSPHSumMassDensity(const ConnectivityMap<Dim<1> >&, 
                                            const TableKernel<Dim<1> >&, 
                                            const FieldList<Dim<1>, Dim<1>::Vector>&,
                                            const FieldList<Dim<1>, Dim<1>::Scalar>&,
                                            const FieldList<Dim<1>, Dim<1>::Scalar>&,
                                            const FieldList<Dim<1>, Dim<1>::SymTensor>&,
                                            const std::vector<BoundarySpace::Boundary<Dim<1> >*>::const_iterator&,
                                            const std::vector<BoundarySpace::Boundary<Dim<1> >*>::const_iterator&,
                                            FieldList<Dim<1>, Dim<1>::Scalar>&);
    template void computeCSPHSumMassDensity(const ConnectivityMap<Dim<2> >&, 
                                            const TableKernel<Dim<2> >&, 
                                            const FieldList<Dim<2>, Dim<2>::Vector>&,
                                            const FieldList<Dim<2>, Dim<2>::Scalar>&,
                                            const FieldList<Dim<2>, Dim<2>::Scalar>&,
                                            const FieldList<Dim<2>, Dim<2>::SymTensor>&,
                                            const std::vector<BoundarySpace::Boundary<Dim<2> >*>::const_iterator&,
                                            const std::vector<BoundarySpace::Boundary<Dim<2> >*>::const_iterator&,
                                            FieldList<Dim<2>, Dim<2>::Scalar>&);
    template void computeCSPHSumMassDensity(const ConnectivityMap<Dim<3> >&, 
                                            const TableKernel<Dim<3> >&, 
                                            const FieldList<Dim<3>, Dim<3>::Vector>&,
                                            const FieldList<Dim<3>, Dim<3>::Scalar>&,
                                            const FieldList<Dim<3>, Dim<3>::Scalar>&,
                                            const FieldList<Dim<3>, Dim<3>::SymTensor>&,
                                            const std::vector<BoundarySpace::Boundary<Dim<3> >*>::const_iterator&,
                                            const std::vector<BoundarySpace::Boundary<Dim<3> >*>::const_iterator&,
                                            FieldList<Dim<3>, Dim<3>::Scalar>&);
  }
}

