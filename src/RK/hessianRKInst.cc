//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "RK/hessianRK.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
template
FieldList<Dim<1>, MathTraits<Dim<1>, Dim<1>::Scalar>::HessianType>
hessianRK<Dim<1>, Dim<1>::Scalar>(const FieldList<Dim<1>, Dim<1>::Scalar>& fieldList,
                                                     const FieldList<Dim<1>, Dim<1>::Vector>& position,
                                                     const FieldList<Dim<1>, Dim<1>::Scalar>& weight,
                                                     const FieldList<Dim<1>, Dim<1>::SymTensor>& H,
                                                     const ConnectivityMap<Dim<1>>& connectivityMap,
                                                     const ReproducingKernel<Dim<1>>& WR,
                                                     const FieldList<Dim<1>, RKCoefficients<Dim<1>>>& corrections,
                                                     const NodeCoupling& nodeCoupling);

template
FieldList<Dim<1>, MathTraits<Dim<1>, Dim<1>::Vector>::HessianType>
hessianRK<Dim<1>, Dim<1>::Vector>(const FieldList<Dim<1>, Dim<1>::Vector>& fieldList,
                                                     const FieldList<Dim<1>, Dim<1>::Vector>& position,
                                                     const FieldList<Dim<1>, Dim<1>::Scalar>& weight,
                                                     const FieldList<Dim<1>, Dim<1>::SymTensor>& H,
                                                     const ConnectivityMap<Dim<1>>& connectivityMap,
                                                     const ReproducingKernel<Dim<1>>& WR,
                                                     const FieldList<Dim<1>, RKCoefficients<Dim<1>>>& corrections,
                                                     const NodeCoupling& nodeCoupling);
#endif

#if defined(SPHERAL_ENABLE_2D)
template
FieldList<Dim<2>, MathTraits<Dim<2>, Dim<2>::Scalar>::HessianType>
hessianRK<Dim<2>, Dim<2>::Scalar>(const FieldList<Dim<2>, Dim<2>::Scalar>& fieldList,
                                                     const FieldList<Dim<2>, Dim<2>::Vector>& position,
                                                     const FieldList<Dim<2>, Dim<2>::Scalar>& weight,
                                                     const FieldList<Dim<2>, Dim<2>::SymTensor>& H,
                                                     const ConnectivityMap<Dim<2>>& connectivityMap,
                                                     const ReproducingKernel<Dim<2>>& WR,
                                                     const FieldList<Dim<2>, RKCoefficients<Dim<2>>>& corrections,
                                                     const NodeCoupling& nodeCoupling);

template
FieldList<Dim<2>, MathTraits<Dim<2>, Dim<2>::Vector>::HessianType>
hessianRK<Dim<2>, Dim<2>::Vector>(const FieldList<Dim<2>, Dim<2>::Vector>& fieldList,
                                                     const FieldList<Dim<2>, Dim<2>::Vector>& position,
                                                     const FieldList<Dim<2>, Dim<2>::Scalar>& weight,
                                                     const FieldList<Dim<2>, Dim<2>::SymTensor>& H,
                                                     const ConnectivityMap<Dim<2>>& connectivityMap,
                                                     const ReproducingKernel<Dim<2>>& WR,
                                                     const FieldList<Dim<2>, RKCoefficients<Dim<2>>>& corrections,
                                                     const NodeCoupling& nodeCoupling);
#endif

#if defined(SPHERAL_ENABLE_3D)
template
FieldList<Dim<3>, MathTraits<Dim<3>, Dim<3>::Scalar>::HessianType>
hessianRK<Dim<3>, Dim<3>::Scalar>(const FieldList<Dim<3>, Dim<3>::Scalar>& fieldList,
                                                     const FieldList<Dim<3>, Dim<3>::Vector>& position,
                                                     const FieldList<Dim<3>, Dim<3>::Scalar>& weight,
                                                     const FieldList<Dim<3>, Dim<3>::SymTensor>& H,
                                                     const ConnectivityMap<Dim<3>>& connectivityMap,
                                                     const ReproducingKernel<Dim<3>>& WR,
                                                     const FieldList<Dim<3>, RKCoefficients<Dim<3>>>& corrections,
                                                     const NodeCoupling& nodeCoupling);

template
FieldList<Dim<3>, MathTraits<Dim<3>, Dim<3>::Vector>::HessianType>
hessianRK<Dim<3>, Dim<3>::Vector>(const FieldList<Dim<3>, Dim<3>::Vector>& fieldList,
                                                     const FieldList<Dim<3>, Dim<3>::Vector>& position,
                                                     const FieldList<Dim<3>, Dim<3>::Scalar>& weight,
                                                     const FieldList<Dim<3>, Dim<3>::SymTensor>& H,
                                                     const ConnectivityMap<Dim<3>>& connectivityMap,
                                                     const ReproducingKernel<Dim<3>>& WR,
                                                     const FieldList<Dim<3>, RKCoefficients<Dim<3>>>& corrections,
                                                     const NodeCoupling& nodeCoupling);
#endif
}
