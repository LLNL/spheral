//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "gradDivVectorFieldList.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  namespace FieldSpace {

    using KernelSpace::TableKernel;

    //============================== gradDivVectorFieldList() ==============================
    template 
    FieldList<Dim<1>, Dim<1>::Vector> 
    gradDivVectorFieldList< Dim<1> >
    (const FieldList<Dim<1>, Dim<1>::Vector>& fieldList,
     const FieldList<Dim<1>, Dim<1>::Vector>& position,
     const FieldList<Dim<1>, Dim<1>::Scalar>& weight,
     const FieldList<Dim<1>, Dim<1>::Scalar>& mass,
     const FieldList<Dim<1>, Dim<1>::Scalar>& density,
     const FieldList<Dim<1>, Dim<1>::SymTensor>& Hfield,
     const TableKernel< Dim<1> >& kernel,
     const vector<Boundary<Dim<1> >*>& boundaries);

    template 
    FieldList<Dim<2>, Dim<2>::Vector> 
    gradDivVectorFieldList< Dim<2> >
    (const FieldList<Dim<2>, Dim<2>::Vector>& fieldList,
     const FieldList<Dim<2>, Dim<2>::Vector>& position,
     const FieldList<Dim<2>, Dim<2>::Scalar>& weight,
     const FieldList<Dim<2>, Dim<2>::Scalar>& mass,
     const FieldList<Dim<2>, Dim<2>::Scalar>& density,
     const FieldList<Dim<2>, Dim<2>::SymTensor>& Hfield,
     const TableKernel< Dim<2> >& kernel,
     const vector<Boundary<Dim<2> >*>& boundaries);

    template 
    FieldList<Dim<3>, Dim<3>::Vector> 
    gradDivVectorFieldList< Dim<3> >
    (const FieldList<Dim<3>, Dim<3>::Vector>& fieldList,
     const FieldList<Dim<3>, Dim<3>::Vector>& position,
     const FieldList<Dim<3>, Dim<3>::Scalar>& weight,
     const FieldList<Dim<3>, Dim<3>::Scalar>& mass,
     const FieldList<Dim<3>, Dim<3>::Scalar>& density,
     const FieldList<Dim<3>, Dim<3>::SymTensor>& Hfield,
     const TableKernel< Dim<3> >& kernel,
     const vector<Boundary<Dim<3> >*>& boundaries);

  }
}
