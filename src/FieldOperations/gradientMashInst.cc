//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "gradientMash.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  namespace FieldSpace {

    using KernelSpace::TableKernel;

    //============================== gradientMash() ==============================
    template 
    FieldList<Dim<1>, MathTraits<Dim<1>, Dim<1>::Scalar>::GradientType> 
    gradientMash<Dim<1>, Dim<1>::Scalar>(const FieldList<Dim<1>, Dim<1>::Scalar>& fieldList,
                                         const FieldList<Dim<1>, Dim<1>::Vector>& position,
                                         const FieldList<Dim<1>, Dim<1>::Scalar>& weight,
                                         const FieldList<Dim<1>, Dim<1>::SymTensor>& Hfield,
                                         const TableKernel< Dim<1> >& kernel);
    template 
    FieldList<Dim<1>, MathTraits<Dim<1>, Dim<1>::Vector>::GradientType> 
    gradientMash<Dim<1>, Dim<1>::Vector>(const FieldList<Dim<1>, Dim<1>::Vector>& fieldList,
                                         const FieldList<Dim<1>, Dim<1>::Vector>& position,
                                         const FieldList<Dim<1>, Dim<1>::Scalar>& weight,
                                         const FieldList<Dim<1>, Dim<1>::SymTensor>& Hfield,
                                         const TableKernel< Dim<1> >& kernel);

    template 
    FieldList<Dim<2>, MathTraits<Dim<2>, Dim<2>::Scalar>::GradientType> 
    gradientMash<Dim<2>, Dim<2>::Scalar>(const FieldList<Dim<2>, Dim<2>::Scalar>& fieldList,
                                         const FieldList<Dim<2>, Dim<2>::Vector>& position,
                                         const FieldList<Dim<2>, Dim<2>::Scalar>& weight,
                                         const FieldList<Dim<2>, Dim<2>::SymTensor>& Hfield,
                                         const TableKernel< Dim<2> >& kernel);
    template 
    FieldList<Dim<2>, MathTraits<Dim<2>, Dim<2>::Vector>::GradientType> 
    gradientMash<Dim<2>, Dim<2>::Vector>(const FieldList<Dim<2>, Dim<2>::Vector>& fieldList,
                                         const FieldList<Dim<2>, Dim<2>::Vector>& position,
                                         const FieldList<Dim<2>, Dim<2>::Scalar>& weight,
                                         const FieldList<Dim<2>, Dim<2>::SymTensor>& Hfield,
                                         const TableKernel< Dim<2> >& kernel);

    template 
    FieldList<Dim<3>, MathTraits<Dim<3>, Dim<3>::Scalar>::GradientType> 
    gradientMash<Dim<3>, Dim<3>::Scalar>(const FieldList<Dim<3>, Dim<3>::Scalar>& fieldList,
                                         const FieldList<Dim<3>, Dim<3>::Vector>& position,
                                         const FieldList<Dim<3>, Dim<3>::Scalar>& weight,
                                         const FieldList<Dim<3>, Dim<3>::SymTensor>& Hfield,
                                         const TableKernel< Dim<3> >& kernel);
    template 
    FieldList<Dim<3>, MathTraits<Dim<3>, Dim<3>::Vector>::GradientType> 
    gradientMash<Dim<3>, Dim<3>::Vector>(const FieldList<Dim<3>, Dim<3>::Vector>& fieldList,
                                         const FieldList<Dim<3>, Dim<3>::Vector>& position,
                                         const FieldList<Dim<3>, Dim<3>::Scalar>& weight,
                                         const FieldList<Dim<3>, Dim<3>::SymTensor>& Hfield,
                                         const TableKernel< Dim<3> >& kernel);


    //============================== gradientMash() ==============================
    template 
    FieldList<Dim<1>, MathTraits<Dim<1>, Dim<1>::Scalar>::GradientType> 
    gradientMash2<Dim<1>, Dim<1>::Scalar>(const FieldList<Dim<1>, Dim<1>::Scalar>& fieldList,
                                          const FieldList<Dim<1>, Dim<1>::Vector>& position,
                                          const FieldList<Dim<1>, Dim<1>::Scalar>& weight,
                                          const FieldList<Dim<1>, Dim<1>::Scalar>& weightDensity,
                                          const FieldList<Dim<1>, Dim<1>::SymTensor>& Hfield,
                                          const TableKernel< Dim<1> >& kernel);

  }
}
