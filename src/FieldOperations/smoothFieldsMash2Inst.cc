//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "smoothFieldsMash2.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  namespace FieldSpace {

    using KernelSpace::TableKernel;

    //============================== smoothFieldsMash2() ==============================
    template 
    FieldList<Dim<1>, Dim<1>::Scalar> 
    smoothFieldsMash2<Dim<1>, Dim<1>::Scalar>(const FieldList<Dim<1>, Dim<1>::Scalar>& fieldList,
                                              const FieldList<Dim<1>, Dim<1>::Vector>& position,
                                              const FieldList<Dim<1>, Dim<1>::Scalar>& weight,
                                              const FieldList<Dim<1>, Dim<1>::Scalar>& weightDensity,
                                              const FieldList<Dim<1>, Dim<1>::SymTensor>& Hfield,
                                              const TableKernel< Dim<1> >& kernel);
    template 
    FieldList<Dim<1>, Dim<1>::Vector> 
    smoothFieldsMash2<Dim<1>, Dim<1>::Vector>(const FieldList<Dim<1>, Dim<1>::Vector>& fieldList,
                                              const FieldList<Dim<1>, Dim<1>::Vector>& position,
                                              const FieldList<Dim<1>, Dim<1>::Scalar>& weight,
                                              const FieldList<Dim<1>, Dim<1>::Scalar>& weightDensity,
                                              const FieldList<Dim<1>, Dim<1>::SymTensor>& Hfield,
                                              const TableKernel< Dim<1> >& kernel);
    template 
    FieldList<Dim<1>, Dim<1>::Tensor> 
    smoothFieldsMash2<Dim<1>, Dim<1>::Tensor>(const FieldList<Dim<1>, Dim<1>::Tensor>& fieldList,
                                              const FieldList<Dim<1>, Dim<1>::Vector>& position,
                                              const FieldList<Dim<1>, Dim<1>::Scalar>& weight,
                                              const FieldList<Dim<1>, Dim<1>::Scalar>& weightDensity,
                                              const FieldList<Dim<1>, Dim<1>::SymTensor>& Hfield,
                                              const TableKernel< Dim<1> >& kernel);
    template 
    FieldList<Dim<1>, Dim<1>::SymTensor> 
    smoothFieldsMash2<Dim<1>, Dim<1>::SymTensor>(const FieldList<Dim<1>, Dim<1>::SymTensor>& fieldList,
                                                 const FieldList<Dim<1>, Dim<1>::Vector>& position,
                                                 const FieldList<Dim<1>, Dim<1>::Scalar>& weight,
                                                 const FieldList<Dim<1>, Dim<1>::Scalar>& weightDensity,
                                                 const FieldList<Dim<1>, Dim<1>::SymTensor>& Hfield,
                                                 const TableKernel< Dim<1> >& kernel);

    template 
    FieldList<Dim<2>, Dim<2>::Scalar> 
    smoothFieldsMash2<Dim<2>, Dim<2>::Scalar>(const FieldList<Dim<2>, Dim<2>::Scalar>& fieldList,
                                              const FieldList<Dim<2>, Dim<2>::Vector>& position,
                                              const FieldList<Dim<2>, Dim<2>::Scalar>& weight,
                                              const FieldList<Dim<2>, Dim<2>::Scalar>& weightDensity,
                                              const FieldList<Dim<2>, Dim<2>::SymTensor>& Hfield,
                                              const TableKernel< Dim<2> >& kernel);
    template 
    FieldList<Dim<2>, Dim<2>::Vector> 
    smoothFieldsMash2<Dim<2>, Dim<2>::Vector>(const FieldList<Dim<2>, Dim<2>::Vector>& fieldList,
                                              const FieldList<Dim<2>, Dim<2>::Vector>& position,
                                              const FieldList<Dim<2>, Dim<2>::Scalar>& weight,
                                              const FieldList<Dim<2>, Dim<2>::Scalar>& weightDensity,
                                              const FieldList<Dim<2>, Dim<2>::SymTensor>& Hfield,
                                              const TableKernel< Dim<2> >& kernel);
    template 
    FieldList<Dim<2>, Dim<2>::Tensor> 
    smoothFieldsMash2<Dim<2>, Dim<2>::Tensor>(const FieldList<Dim<2>, Dim<2>::Tensor>& fieldList,
                                              const FieldList<Dim<2>, Dim<2>::Vector>& position,
                                              const FieldList<Dim<2>, Dim<2>::Scalar>& weight,
                                              const FieldList<Dim<2>, Dim<2>::Scalar>& weightDensity,
                                              const FieldList<Dim<2>, Dim<2>::SymTensor>& Hfield,
                                              const TableKernel< Dim<2> >& kernel);
    template 
    FieldList<Dim<2>, Dim<2>::SymTensor> 
    smoothFieldsMash2<Dim<2>, Dim<2>::SymTensor>(const FieldList<Dim<2>, Dim<2>::SymTensor>& fieldList,
                                                 const FieldList<Dim<2>, Dim<2>::Vector>& position,
                                                 const FieldList<Dim<2>, Dim<2>::Scalar>& weight,
                                                 const FieldList<Dim<2>, Dim<2>::Scalar>& weightDensity,
                                                 const FieldList<Dim<2>, Dim<2>::SymTensor>& Hfield,
                                                 const TableKernel< Dim<2> >& kernel);

    template 
    FieldList<Dim<3>, Dim<3>::Scalar> 
    smoothFieldsMash2<Dim<3>, Dim<3>::Scalar>(const FieldList<Dim<3>, Dim<3>::Scalar>& fieldList,
                                              const FieldList<Dim<3>, Dim<3>::Vector>& position,
                                              const FieldList<Dim<3>, Dim<3>::Scalar>& weight,
                                              const FieldList<Dim<3>, Dim<3>::Scalar>& weightDensity,
                                              const FieldList<Dim<3>, Dim<3>::SymTensor>& Hfield,
                                              const TableKernel< Dim<3> >& kernel);
    template 
    FieldList<Dim<3>, Dim<3>::Vector> 
    smoothFieldsMash2<Dim<3>, Dim<3>::Vector>(const FieldList<Dim<3>, Dim<3>::Vector>& fieldList,
                                              const FieldList<Dim<3>, Dim<3>::Vector>& position,
                                              const FieldList<Dim<3>, Dim<3>::Scalar>& weight,
                                              const FieldList<Dim<3>, Dim<3>::Scalar>& weightDensity,
                                              const FieldList<Dim<3>, Dim<3>::SymTensor>& Hfield,
                                              const TableKernel< Dim<3> >& kernel);
    template 
    FieldList<Dim<3>, Dim<3>::Tensor> 
    smoothFieldsMash2<Dim<3>, Dim<3>::Tensor>(const FieldList<Dim<3>, Dim<3>::Tensor>& fieldList,
                                              const FieldList<Dim<3>, Dim<3>::Vector>& position,
                                              const FieldList<Dim<3>, Dim<3>::Scalar>& weight,
                                              const FieldList<Dim<3>, Dim<3>::Scalar>& weightDensity,
                                              const FieldList<Dim<3>, Dim<3>::SymTensor>& Hfield,
                                              const TableKernel< Dim<3> >& kernel);
    template 
    FieldList<Dim<3>, Dim<3>::SymTensor> 
    smoothFieldsMash2<Dim<3>, Dim<3>::SymTensor>(const FieldList<Dim<3>, Dim<3>::SymTensor>& fieldList,
                                                 const FieldList<Dim<3>, Dim<3>::Vector>& position,
                                                 const FieldList<Dim<3>, Dim<3>::Scalar>& weight,
                                                 const FieldList<Dim<3>, Dim<3>::Scalar>& weightDensity,
                                                 const FieldList<Dim<3>, Dim<3>::SymTensor>& Hfield,
                                                 const TableKernel< Dim<3> >& kernel);
  }
}
