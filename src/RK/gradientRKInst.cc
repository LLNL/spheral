//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "RK/gradientRK.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
"""

for Value in ("Scalar", "Vector"):
    text += """
template
FieldList<Dim<%1>, MathTraits<Dim<%1>, Dim<%1>::%(Value)s>::GradientType>
gradientRK<Dim<%1>, Dim<%1>::%(Value)s>(const FieldList<Dim<%1>, Dim<%1>::%(Value)s>& fieldList,
                                                      const FieldList<Dim<%1>, Dim<%1>::Vector>& position,
                                                      const FieldList<Dim<%1>, Dim<%1>::Scalar>& weight,
                                                      const FieldList<Dim<%1>, Dim<%1>::SymTensor>& H,
                                                      const ConnectivityMap<Dim<%1>>& connectivityMap,
                                                      const ReproducingKernel<Dim<%1>>& WR,
                                                      const FieldList<Dim<%1>, RKCoefficients<Dim<%1>>>& corrections,
                                                      const NodeCoupling& nodeCoupling);

""" % {"Value" : Value}

    text += """
template
FieldList<Dim<%1>, std::vector<MathTraits<Dim<%1>, Dim<%1>::%(Value)s>::GradientType>>
gradientRK<Dim<%1>, Dim<%1>::%(Value)s>(const FieldList<Dim<%1>, std::vector<Dim<%1>::%(Value)s>>& fieldList,
                                                      const FieldList<Dim<%1>, Dim<%1>::Vector>& position,
                                                      const FieldList<Dim<%1>, Dim<%1>::Scalar>& weight,
                                                      const FieldList<Dim<%1>, Dim<%1>::SymTensor>& H,
                                                      const ConnectivityMap<Dim<%1>>& connectivityMap,
                                                      const ReproducingKernel<Dim<%1>>& WR,
                                                      const FieldList<Dim<%1>, RKCoefficients<Dim<%1>>>& corrections,
                                                      const NodeCoupling& nodeCoupling);

""" % {"Value" : Value}

text += """
#endif

#if defined(SPHERAL_ENABLE_2D)
"""

for Value in ("Scalar", "Vector"):
    text += """
template
FieldList<Dim<%2>, MathTraits<Dim<%2>, Dim<%2>::%(Value)s>::GradientType>
gradientRK<Dim<%2>, Dim<%2>::%(Value)s>(const FieldList<Dim<%2>, Dim<%2>::%(Value)s>& fieldList,
                                                      const FieldList<Dim<%2>, Dim<%2>::Vector>& position,
                                                      const FieldList<Dim<%2>, Dim<%2>::Scalar>& weight,
                                                      const FieldList<Dim<%2>, Dim<%2>::SymTensor>& H,
                                                      const ConnectivityMap<Dim<%2>>& connectivityMap,
                                                      const ReproducingKernel<Dim<%2>>& WR,
                                                      const FieldList<Dim<%2>, RKCoefficients<Dim<%2>>>& corrections,
                                                      const NodeCoupling& nodeCoupling);

""" % {"Value" : Value}

    text += """
template
FieldList<Dim<%2>, std::vector<MathTraits<Dim<%2>, Dim<%2>::%(Value)s>::GradientType>>
gradientRK<Dim<%2>, Dim<%2>::%(Value)s>(const FieldList<Dim<%2>, std::vector<Dim<%2>::%(Value)s>>& fieldList,
                                                      const FieldList<Dim<%2>, Dim<%2>::Vector>& position,
                                                      const FieldList<Dim<%2>, Dim<%2>::Scalar>& weight,
                                                      const FieldList<Dim<%2>, Dim<%2>::SymTensor>& H,
                                                      const ConnectivityMap<Dim<%2>>& connectivityMap,
                                                      const ReproducingKernel<Dim<%2>>& WR,
                                                      const FieldList<Dim<%2>, RKCoefficients<Dim<%2>>>& corrections,
                                                      const NodeCoupling& nodeCoupling);

""" % {"Value" : Value}

text += """
#endif

#if defined(SPHERAL_ENABLE_3D)
"""

for Value in ("Scalar", "Vector"):
    text += """
template
FieldList<Dim<%3>, MathTraits<Dim<%3>, Dim<%3>::%(Value)s>::GradientType>
gradientRK<Dim<%3>, Dim<%3>::%(Value)s>(const FieldList<Dim<%3>, Dim<%3>::%(Value)s>& fieldList,
                                                      const FieldList<Dim<%3>, Dim<%3>::Vector>& position,
                                                      const FieldList<Dim<%3>, Dim<%3>::Scalar>& weight,
                                                      const FieldList<Dim<%3>, Dim<%3>::SymTensor>& H,
                                                      const ConnectivityMap<Dim<%3>>& connectivityMap,
                                                      const ReproducingKernel<Dim<%3>>& WR,
                                                      const FieldList<Dim<%3>, RKCoefficients<Dim<%3>>>& corrections,
                                                      const NodeCoupling& nodeCoupling);

""" % {"Value" : Value}

    text += """
template
FieldList<Dim<%3>, std::vector<MathTraits<Dim<%3>, Dim<%3>::%(Value)s>::GradientType>>
gradientRK<Dim<%3>, Dim<%3>::%(Value)s>(const FieldList<Dim<%3>, std::vector<Dim<%3>::%(Value)s>>& fieldList,
                                                      const FieldList<Dim<%3>, Dim<%3>::Vector>& position,
                                                      const FieldList<Dim<%3>, Dim<%3>::Scalar>& weight,
                                                      const FieldList<Dim<%3>, Dim<%3>::SymTensor>& H,
                                                      const ConnectivityMap<Dim<%3>>& connectivityMap,
                                                      const ReproducingKernel<Dim<%3>>& WR,
                                                      const FieldList<Dim<%3>, RKCoefficients<Dim<%3>>>& corrections,
                                                      const NodeCoupling& nodeCoupling);

""" % {"Value" : Value}

text += """
#endif
}