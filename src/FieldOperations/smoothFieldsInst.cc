//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "smoothFields.cc"
#include "Geometry/Dimension.hh"

//============================== smoothFields() ==============================
namespace Spheral {
  namespace FieldSpace {

    using KernelSpace::TableKernel;

    template 
    FieldList<Dim<1>, Dim<1>::Scalar> 
    smoothFields<Dim<1>, Dim<1>::Scalar>(const FieldList<Dim<1>, Dim<1>::Scalar>& fieldList,
                                         const FieldList<Dim<1>, Dim<1>::Vector>& position,
                                         const FieldList<Dim<1>, Dim<1>::Scalar>& weight,
                                         const FieldList<Dim<1>, Dim<1>::Scalar>& mass,
                                         const FieldList<Dim<1>, Dim<1>::Scalar>& density,
                                         const FieldList<Dim<1>, Dim<1>::SymTensor>& Hfield,
                                         const TableKernel< Dim<1> >& kernel);
    template 
    FieldList<Dim<1>, Dim<1>::Vector> 
    smoothFields<Dim<1>, Dim<1>::Vector>(const FieldList<Dim<1>, Dim<1>::Vector>& fieldList,
                                         const FieldList<Dim<1>, Dim<1>::Vector>& position,
                                         const FieldList<Dim<1>, Dim<1>::Scalar>& weight,
                                         const FieldList<Dim<1>, Dim<1>::Scalar>& mass,
                                         const FieldList<Dim<1>, Dim<1>::Scalar>& density,
                                         const FieldList<Dim<1>, Dim<1>::SymTensor>& Hfield,
                                         const TableKernel< Dim<1> >& kernel);
    template 
    FieldList<Dim<1>, Dim<1>::Tensor> 
    smoothFields<Dim<1>, Dim<1>::Tensor>(const FieldList<Dim<1>, Dim<1>::Tensor>& fieldList,
                                         const FieldList<Dim<1>, Dim<1>::Vector>& position,
                                         const FieldList<Dim<1>, Dim<1>::Scalar>& weight,
                                         const FieldList<Dim<1>, Dim<1>::Scalar>& mass,
                                         const FieldList<Dim<1>, Dim<1>::Scalar>& density,
                                         const FieldList<Dim<1>, Dim<1>::SymTensor>& Hfield,
                                         const TableKernel< Dim<1> >& kernel);
    template 
    FieldList<Dim<1>, Dim<1>::SymTensor> 
    smoothFields<Dim<1>, Dim<1>::SymTensor>(const FieldList<Dim<1>, Dim<1>::SymTensor>& fieldList,
                                            const FieldList<Dim<1>, Dim<1>::Vector>& position,
                                            const FieldList<Dim<1>, Dim<1>::Scalar>& weight,
                                            const FieldList<Dim<1>, Dim<1>::Scalar>& mass,
                                            const FieldList<Dim<1>, Dim<1>::Scalar>& density,
                                            const FieldList<Dim<1>, Dim<1>::SymTensor>& Hfield,
                                            const TableKernel< Dim<1> >& kernel);

    template 
    FieldList<Dim<2>, Dim<2>::Scalar> 
    smoothFields<Dim<2>, Dim<2>::Scalar>(const FieldList<Dim<2>, Dim<2>::Scalar>& fieldList,
                                         const FieldList<Dim<2>, Dim<2>::Vector>& position,
                                         const FieldList<Dim<2>, Dim<2>::Scalar>& weight,
                                         const FieldList<Dim<2>, Dim<2>::Scalar>& mass,
                                         const FieldList<Dim<2>, Dim<2>::Scalar>& density,
                                         const FieldList<Dim<2>, Dim<2>::SymTensor>& Hfield,
                                         const TableKernel< Dim<2> >& kernel);
    template 
    FieldList<Dim<2>, Dim<2>::Vector> 
    smoothFields<Dim<2>, Dim<2>::Vector>(const FieldList<Dim<2>, Dim<2>::Vector>& fieldList,
                                         const FieldList<Dim<2>, Dim<2>::Vector>& position,
                                         const FieldList<Dim<2>, Dim<2>::Scalar>& weight,
                                         const FieldList<Dim<2>, Dim<2>::Scalar>& mass,
                                         const FieldList<Dim<2>, Dim<2>::Scalar>& density,
                                         const FieldList<Dim<2>, Dim<2>::SymTensor>& Hfield,
                                         const TableKernel< Dim<2> >& kernel);
    template 
    FieldList<Dim<2>, Dim<2>::Tensor> 
    smoothFields<Dim<2>, Dim<2>::Tensor>(const FieldList<Dim<2>, Dim<2>::Tensor>& fieldList,
                                         const FieldList<Dim<2>, Dim<2>::Vector>& position,
                                         const FieldList<Dim<2>, Dim<2>::Scalar>& weight,
                                         const FieldList<Dim<2>, Dim<2>::Scalar>& mass,
                                         const FieldList<Dim<2>, Dim<2>::Scalar>& density,
                                         const FieldList<Dim<2>, Dim<2>::SymTensor>& Hfield,
                                         const TableKernel< Dim<2> >& kernel);
    template 
    FieldList<Dim<2>, Dim<2>::SymTensor> 
    smoothFields<Dim<2>, Dim<2>::SymTensor>(const FieldList<Dim<2>, Dim<2>::SymTensor>& fieldList,
                                            const FieldList<Dim<2>, Dim<2>::Vector>& position,
                                            const FieldList<Dim<2>, Dim<2>::Scalar>& weight,
                                            const FieldList<Dim<2>, Dim<2>::Scalar>& mass,
                                            const FieldList<Dim<2>, Dim<2>::Scalar>& density,
                                            const FieldList<Dim<2>, Dim<2>::SymTensor>& Hfield,
                                            const TableKernel< Dim<2> >& kernel);

    template 
    FieldList<Dim<3>, Dim<3>::Scalar> 
    smoothFields<Dim<3>, Dim<3>::Scalar>(const FieldList<Dim<3>, Dim<3>::Scalar>& fieldList,
                                         const FieldList<Dim<3>, Dim<3>::Vector>& position,
                                         const FieldList<Dim<3>, Dim<3>::Scalar>& weight,
                                         const FieldList<Dim<3>, Dim<3>::Scalar>& mass,
                                         const FieldList<Dim<3>, Dim<3>::Scalar>& density,
                                         const FieldList<Dim<3>, Dim<3>::SymTensor>& Hfield,
                                         const TableKernel< Dim<3> >& kernel);
    template 
    FieldList<Dim<3>, Dim<3>::Vector> 
    smoothFields<Dim<3>, Dim<3>::Vector>(const FieldList<Dim<3>, Dim<3>::Vector>& fieldList,
                                         const FieldList<Dim<3>, Dim<3>::Vector>& position,
                                         const FieldList<Dim<3>, Dim<3>::Scalar>& weight,
                                         const FieldList<Dim<3>, Dim<3>::Scalar>& mass,
                                         const FieldList<Dim<3>, Dim<3>::Scalar>& density,
                                         const FieldList<Dim<3>, Dim<3>::SymTensor>& Hfield,
                                         const TableKernel< Dim<3> >& kernel);
    template 
    FieldList<Dim<3>, Dim<3>::Tensor> 
    smoothFields<Dim<3>, Dim<3>::Tensor>(const FieldList<Dim<3>, Dim<3>::Tensor>& fieldList,
                                         const FieldList<Dim<3>, Dim<3>::Vector>& position,
                                         const FieldList<Dim<3>, Dim<3>::Scalar>& weight,
                                         const FieldList<Dim<3>, Dim<3>::Scalar>& mass,
                                         const FieldList<Dim<3>, Dim<3>::Scalar>& density,
                                         const FieldList<Dim<3>, Dim<3>::SymTensor>& Hfield,
                                         const TableKernel< Dim<3> >& kernel);
    template 
    FieldList<Dim<3>, Dim<3>::SymTensor> 
    smoothFields<Dim<3>, Dim<3>::SymTensor>(const FieldList<Dim<3>, Dim<3>::SymTensor>& fieldList,
                                            const FieldList<Dim<3>, Dim<3>::Vector>& position,
                                            const FieldList<Dim<3>, Dim<3>::Scalar>& weight,
                                            const FieldList<Dim<3>, Dim<3>::Scalar>& mass,
                                            const FieldList<Dim<3>, Dim<3>::Scalar>& density,
                                            const FieldList<Dim<3>, Dim<3>::SymTensor>& Hfield,
                                            const TableKernel< Dim<3> >& kernel);
  }
}
