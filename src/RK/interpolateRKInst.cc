//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "RK/interpolateRK.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)

template
std::vector<boost::variant<FieldList<Dim<1>, Dim<1>::Scalar>,
                           FieldList<Dim<1>, Dim<1>::Vector>,
                           FieldList<Dim<1>, Dim<1>::Tensor>,
                           FieldList<Dim<1>, Dim<1>::SymTensor>,
                           FieldList<Dim<1>, Dim<1>::ThirdRankTensor>>>
interpolateRK<Dim<1>>(const std::vector<boost::variant<FieldList<Dim<1>, Dim<1>::Scalar>,
                                                              FieldList<Dim<1>, Dim<1>::Vector>,
                                                              FieldList<Dim<1>, Dim<1>::Tensor>,
                                                              FieldList<Dim<1>, Dim<1>::SymTensor>,
                                                              FieldList<Dim<1>, Dim<1>::ThirdRankTensor>>>& fieldLists,
                                          const FieldList<Dim<1>, Dim<1>::Vector>& position,
                                          const FieldList<Dim<1>, Dim<1>::Scalar>& weight,
                                          const FieldList<Dim<1>, Dim<1>::SymTensor>& H,
                                          const ConnectivityMap<Dim<1> >& connectivityMap,
                                          const ReproducingKernel< Dim<1> >& kernel,
                                          const FieldList<Dim<1>, RKCoefficients<Dim<1>>>& corrections,
                                          const NodeCoupling& nodeCoupling);

#endif

#if defined(SPHERAL_ENABLE_2D)

template
std::vector<boost::variant<FieldList<Dim<2>, Dim<2>::Scalar>,
                           FieldList<Dim<2>, Dim<2>::Vector>,
                           FieldList<Dim<2>, Dim<2>::Tensor>,
                           FieldList<Dim<2>, Dim<2>::SymTensor>,
                           FieldList<Dim<2>, Dim<2>::ThirdRankTensor>>>
interpolateRK<Dim<2>>(const std::vector<boost::variant<FieldList<Dim<2>, Dim<2>::Scalar>,
                                                              FieldList<Dim<2>, Dim<2>::Vector>,
                                                              FieldList<Dim<2>, Dim<2>::Tensor>,
                                                              FieldList<Dim<2>, Dim<2>::SymTensor>,
                                                              FieldList<Dim<2>, Dim<2>::ThirdRankTensor>>>& fieldLists,
                                          const FieldList<Dim<2>, Dim<2>::Vector>& position,
                                          const FieldList<Dim<2>, Dim<2>::Scalar>& weight,
                                          const FieldList<Dim<2>, Dim<2>::SymTensor>& H,
                                          const ConnectivityMap<Dim<2> >& connectivityMap,
                                          const ReproducingKernel< Dim<2> >& kernel,
                                          const FieldList<Dim<2>, RKCoefficients<Dim<2>>>& corrections,
                                          const NodeCoupling& nodeCoupling);

#endif

#if defined(SPHERAL_ENABLE_3D)

template
std::vector<boost::variant<FieldList<Dim<3>, Dim<3>::Scalar>,
                           FieldList<Dim<3>, Dim<3>::Vector>,
                           FieldList<Dim<3>, Dim<3>::Tensor>,
                           FieldList<Dim<3>, Dim<3>::SymTensor>,
                           FieldList<Dim<3>, Dim<3>::ThirdRankTensor>>>
interpolateRK<Dim<3>>(const std::vector<boost::variant<FieldList<Dim<3>, Dim<3>::Scalar>,
                                                              FieldList<Dim<3>, Dim<3>::Vector>,
                                                              FieldList<Dim<3>, Dim<3>::Tensor>,
                                                              FieldList<Dim<3>, Dim<3>::SymTensor>,
                                                              FieldList<Dim<3>, Dim<3>::ThirdRankTensor>>>& fieldLists,
                                          const FieldList<Dim<3>, Dim<3>::Vector>& position,
                                          const FieldList<Dim<3>, Dim<3>::Scalar>& weight,
                                          const FieldList<Dim<3>, Dim<3>::SymTensor>& H,
                                          const ConnectivityMap<Dim<3> >& connectivityMap,
                                          const ReproducingKernel< Dim<3> >& kernel,
                                          const FieldList<Dim<3>, RKCoefficients<Dim<3>>>& corrections,
                                          const NodeCoupling& nodeCoupling);

#endif
}